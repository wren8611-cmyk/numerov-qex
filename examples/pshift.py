import numpy as np
import matplotlib.pyplot as plt

import configuration as conf
import numerov_qex as nv

psys = conf.system

channel0 = psys.channel_by_id(0)
channel1 = psys.channel_by_id(1)

print(channel0.mu)
print(channel1.mu)

energies = np.linspace(2.00000001, 4.0, 500)

solver = nv.NumerovSolver(
    conf.system,
    conf.potential,
    rmax=20.0,
    npts=2048,
    bc=nv.BoundaryCondition.Outgoing,
)

S00 = []
S11 = []

for E in energies:
    result = solver.solve(E)

    rsys = result.system
    print(rsys.channel_by_id(0).mu)
    print(rsys.channel_by_id(1).mu)

    S00.append(result.s_flux(0, 0))
    S11.append(result.s_flux(1, 1))

S00 = np.asarray(S00)
S11 = np.asarray(S11)

# Valid open-channel regions.
mask0 = np.isfinite(S00)
mask1 = np.isfinite(S11)

# Phase shifts.
delta0 = np.full(len(energies), np.nan)
delta1 = np.full(len(energies), np.nan)

delta0[mask0] = 0.5 * np.unwrap(np.angle(S00[mask0]))
delta1[mask1] = 0.5 * np.unwrap(np.angle(S11[mask1]))

# Inelasticity.
eta0 = np.abs(S00)
eta1 = np.abs(S11)

# Plot phase shifts.
fig, ax = plt.subplots()

ax.plot(energies, np.degrees(delta0), label=r"$\delta_0$")
ax.plot(energies, np.degrees(delta1), label=r"$\delta_1$")

ax.set_xlabel("Energy")
ax.set_ylabel("Phase shift [deg]")
ax.legend()
ax.grid()

# Plot inelasticity.
fig, ax = plt.subplots()

ax.plot(energies, eta0, label=r"$|S_{00}|$")
ax.plot(energies, eta1, "--", label=r"$|S_{11}|$")

ax.set_xlabel("Energy")
ax.set_ylabel(r"Inelasticity $\eta$")
ax.set_ylim(0, 1.05)
ax.legend()
ax.grid()

plt.show()