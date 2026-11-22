# numerov-qex

**numerov-qex** is a Numerov solver for the radial Schrödinger equation in three dimensions. It computes single- and coupled-channel scattering solutions on a uniform radial grid and returns the radial wavefunctions and S-matrix. The solver is written in C++17 and exposed to Python through pybind11.

## Dependencies

- Python 3.10 or later, with development headers for the selected interpreter.
- NumPy and SciPy. The correlation helper uses `numpy.trapezoid`, which requires NumPy 2.0 or later.
- A C++17 compiler, Make, Eigen, and the GNU Scientific Library (GSL).
- The build backend installs setuptools ≥ 68, wheel, and pybind11 ≥ 2.13 through pip's default build isolation.

## Installation

Clone the repository and create a machine-specific build configuration.

```sh
git clone https://github.com/wren8611-cmyk/numerov-qex.git
cd numerov-qex
cp config.mk.example config.mk
```

In `config.mk`, set `CXX` to your C++17 compiler, `EIGEN_DIR` to the directory containing `Eigen/`, and `GSL_DIR` to the GSL installation prefix containing `include/` and `lib/`.

For example, on macOS with Homebrew, install the libraries with `brew install eigen gsl`. Use the output of `brew --prefix eigen` followed by `/include/eigen3` for `EIGEN_DIR`, and the output of `brew --prefix gsl` for `GSL_DIR`.

From your Python environment, install the package and run the tutorial.

```sh
python -m pip install -e .
python examples/numerov_qex_tutorial.py
```

Use `python -m pip install .` for a regular installation. After changing the C++ compiler or library paths, run `make clean` before reinstalling to force a rebuild.

## Usage

Use consistent natural units with $\hbar=c=1$. Masses, center-of-mass energies, momenta, potential, and Radii all use common units.

### Define the channels and potential

Each channel specifies a unique ID, a name, the two particle masses, orbital angular momentum `ell`, and charges `z1` and `z2` in units of the elementary charge. The threshold is the sum of the two masses.

```python
import numerov_qex as nv
from numerov_qex.observables import inelasticity, phase_shift

ch0, ch1 = 10, 20
# nv.ChannelInfo(channel_id, name_str, m1, m2, ell, z1, z2)
system = nv.ChannelSystem([
    nv.ChannelInfo(ch0, "channel 0", 1.0, 1.0, ell=0, z1=0.0, z2=0.0),
    nv.ChannelInfo(ch1, "channel 1", 1.0, 1.5, ell=0, z1=0.0, z2=0.0),
])

potential = nv.PotentialMatrix([], system)
potential.append([
    nv.PotentialInfo("gaussian", [-5.0, 1.0], ch0, ch0),
    nv.PotentialInfo("gaussian", [-3.0, 1.0], ch1, ch1),
    nv.PotentialInfo("gaussian", [-1.0, 1.0], ch0, ch1),
    nv.PotentialInfo("gaussian", [-1.0, 1.0], ch1, ch0),
])
```

| Potential | Parameters | Radial form |
| --- | --- | --- |
| `"gaussian"` | `[a, b]` | $a~\exp[-(r/b)^2]$ |
| `"well"` | `[a, b]` | $a$ for $r<b$, zero otherwise |

Each term is assigned to a row and column by channel ID. For Hermitian coupling, supply both off-diagonal entries. The solver does not check Hermiticity. The solver adds diagonal Coulomb interactions automatically when the charge product is nonzero. `potential.evaluate(r)` returns the explicit matrix of user-supplied terms at radius `r`.

`append()` accepts one term or a list of terms and adds them to the current potential. Use `clear()` to reset those terms. 
An existing solver sees changes to the same `PotentialMatrix` and uses the updated potential when calling `solve()`.

### Configure the solver and solve

```python
solver = nv.NumerovSolver(
    system,
    potential,
    rmax=20.0,
    npts=2048,
    bc=nv.BoundaryCondition.Outgoing,
)

result = solver.solve(2.7)

# potential.clear() ### for resetting potential
```

The grid includes both endpoints, from zero to `rmax`. Choose `rmax` beyond the range of the short-range potential and check convergence in `rmax` and `npts`. Use outgoing boundary conditions for the phase shift and inelasticity helpers and incoming boundary conditions for the correlation function helper.

Pass the total center-of-mass energy, including rest masses, to `solve(E)`. Both channels are open at `E=2.7`. At `E=2.2`, channel 20 is closed. The matching routine requires energies away from exact channel thresholds. Reuse the solver to scan other energies.

### Read the scattering result

| Attribute or method | What it gives you |
| --- | --- |
| `result.E`, `result.boundary_condition`, `result.system` | Energy, matching convention, and channel system |
| `result.open_channel_ids`, `result.is_open(id)` | Open-channel IDs in S-matrix order and a membership check |
| `result.S`, `result.S_flux` | S-matrix and flux-normalized S-matrix over the open channels |
| `result.open_channel_index(id)` | Array index of an open channel |
| `result.s(row_id, col_id)`, `result.s_flux(row_id, col_id)` | S-matrix elements addressed by channel ID |
| `result.wavefunction.x()`, `result.wavefunction[i]` | Radial grid and the wavefunction matrix at grid index `i` |
| `result.psi(i, component_id, reference_id)` | Reduced radial wavefunction for one component and one open reference channel. The component channel may be open or closed. |

```python
print(result.open_channel_ids)
print(result.S_flux)
print("Phase shift [rad]:", phase_shift(result, ch0))
print("Inelasticity:", inelasticity(result, ch0))

radii = result.wavefunction.x()
u_ch1_from_ch0 = [result.psi(i, ch1, ch0) for i in range(len(radii))]
```

The helpers return $\delta_i=\arg S_{ii}/2$ and $\eta_i=|S_{ii}|$. 
Phase shifts are in radians on the principal branch $(-\pi,\pi)$. Both helpers return `nan` for a closed channel.

### Compute a femtoscopic correlation contribution

Construct an incoming-boundary solver and supply a normalized radial source $S(r)$ with $4\pi\int_0^\infty r^2 S(r)dr=1$. The channel IDs in the weights dictionary must have the same orbital angular momentum as the incoming channel.

```python
import numpy as np
from numerov_qex.observables import correlation_function

incoming_solver = nv.NumerovSolver(
    system, potential, bc=nv.BoundaryCondition.Incoming
)

R = 1.0

def source(r):
    return np.exp(-r**2 / (4.0 * R**2)) / (4.0 * np.pi * R**2)**1.5

q, delta_c0 = correlation_function(
    incoming_solver.solve(2.7),
    source,
    weights={ch0: 1.0, ch1: 1.0},
    incoming_channel=ch0,
    ell=0,
)
c_s_wave = 1.0 + delta_c0
```

The helper returns the incoming-channel relative momentum $q$ and the partial-wave correction $\Delta C_\ell(q)$. Combine partial waves at a common $q$ using $C(q)=1+\sum_\ell(2\ell+1)\Delta C_\ell(q)$.

See [the tutorial](examples/numerov_qex_tutorial.py) for a runnable walkthrough of the scattering API.

## AI disclosure

This solver was developed with the assistance of ChatGPT (GPT-6). 
Simple tests (free, free+coulomb, well in single- and coupled-channel) have been performed by a human.

## Version

**0.1.0** is the initial version of numerov-qex. The package version is defined in [`pyproject.toml`](pyproject.toml).
