import numpy as np
import matplotlib.pyplot as plt

import configuration as conf
import numerov_qex as nv

from numerov_qex.observables import correlation_function


# Solver
solver = nv.NumerovSolver(
    conf.system,
    conf.potential,
    rmax=20.0,
    npts=2048,
    bc=nv.BoundaryCondition.Incoming,
)

# Parameters
R = 1.0
incoming_channel = 0
ell = 0

weights = {
    0: 1.0,
    1: 1.0,
}

# Gaussian source
def source(r):
    return np.exp(-r**2 / (4.0 * R**2)) / (4.0 * np.pi * R**2)**1.5


# Incoming channel information
channel = conf.system.channel_by_id(incoming_channel)
mu = channel.mu
threshold = channel.threshold

# Momentum grid
qgrid = np.linspace(0.01, 3.0, 500)

# Correlation function
q_values = []
C_values = []

for q in qgrid:
    E = threshold + q**2 / (2.0 * mu)

    result = solver.solve(E)

    q_out, C = correlation_function(
        result,
        source=source,
        weights=weights,
        incoming_channel=incoming_channel,
        ell=ell,
    )

    q_values.append(q_out)
    C_values.append(C)

q_values = np.asarray(q_values)
C_values = np.asarray(C_values)

# Plot
plt.figure(figsize=(8, 5))

plt.plot(q_values, C_values, color="crimson", zorder=99)

plt.axhline(1.0, color="black", linestyle="--")
plt.xlabel("q")
plt.ylabel(rf"$C_{{{ell}}}(q)$")
plt.grid()
plt.show()