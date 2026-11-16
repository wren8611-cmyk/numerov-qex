# femto.py

import numpy as np
from scipy.special import spherical_jn
import matplotlib.pyplot as plt

import configuration as conf
import numerov_qex as nv

solver = nv.NumerovSolver(
    conf.system,
    conf.potential,
    rmax=20.0,
    npts=2048,
    bc=nv.BoundaryCondition.Incoming,
)

# Parameters
R = 1.0
w = [1.0, 1.0]

# Gaussian source
def source(r, R):
    return np.exp(-r**2 / (4.0 * R**2)) / (4.0 * np.pi * R**2)**1.5

# Observed channel 0
mu0 = 0.5
threshold0 = 2.0

# Momentum grid
qgrid = np.linspace(0.01, 3.0, 500)

C = []

for q in qgrid:

    E = threshold0 + q**2 / (2.0 * mu0)

    result = solver.solve(E)

    # Normalize using incoming channel 0
    norm = 2j * q

    # Avoid r=0
    r = solver.rgrid[1:]
    
    correction = np.zeros_like(r)
    
    for ch in range(2):
    
        u = np.array([
            result.psi(i, ch, 0)
            for i in range(1, solver.npts)
        ])
    
        psi = u / (norm * r)
    
        correction += w[ch] * np.abs(psi)**2

    # Subtract free S-wave
    correction -= spherical_jn(0, q * r)**2

    Cq = 1.0 + 4.0 * np.pi * np.trapezoid(
        r**2 * source(r, R) * correction,
        r
    )

    C.append(Cq)

C = np.asarray(C)

plt.figure(figsize=(8, 5))
plt.plot(qgrid, C, zorder=99, color="crimson")

plt.axhline(1.0, color="black", linestyle="--")
plt.xlabel("q")
plt.ylabel("C(q)")
plt.grid()
plt.show()