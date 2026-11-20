### Dynamic Potential Updates

`NumerovSolver` retains a reference to the `PotentialMatrix` supplied during construction. Therefore, modifications to the potential are automatically reflected in subsequent calls to `solve()`, without reconstructing the solver.

```python
import numerov_qex.radialsolver as rs
import numerov_qex.solver as nv

system = rs.ChannelSystem([
    rs.ChannelInfo(0, "channel0", 1.0, 1.0),
    rs.ChannelInfo(1, "channel1", 1.0, 1.5),
])

potential = rs.PotentialMatrix([], system)

potential.append([
    rs.PotentialInfo("gaussian", [-5.0, 1.0], 0, 0),
    rs.PotentialInfo("gaussian", [-3.0, 1.0], 1, 1),
])

solver = nv.NumerovSolver(system=system, potential=potential)

# Solve using the initial potential.
result1 = solver.solve(3.0)

# Modify the existing potential.
potential.append(
    rs.PotentialInfo("gaussian", [-2.0, 1.0], 0, 0)
)

# The solver automatically uses the updated potential.
result2 = solver.solve(3.0)
```

The two results correspond to different potentials, even though the same solver instance is used.

The potential can also be cleared and reconstructed:

```python
potential.clear()
potential.append(new_potentials)

result3 = solver.solve(3.0)
```

This behavior is useful for parameter scans and fitting procedures, where the potential changes repeatedly while the channel configuration remains fixed.

**Important notes:**

- Modifications to `PotentialMatrix` affect subsequent calls to `solve()`.
- The solver does not need to be reconstructed after modifying the potential.
- `clear()` removes all potential terms without changing the channel configuration.
- Modifying the potential during an active `solve()` call is not supported.
- Previously computed results correspond to the potential used at the time of their calculation.