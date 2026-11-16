# solver.py
# c++ numerov solver interface

import numpy as np
from . import radialsolver as rs

class NumerovSolver:
    def __init__(
        self,
        system,
        potential,
        *,
        rmax=20.0,
        npts=2048,
        bc=rs.BoundaryCondition.Outgoing,
    ):
        if not np.isfinite(rmax) or rmax <= 0:
            raise ValueError("rmax must be finite and positive")

        if not isinstance(npts, (int, np.integer)) or isinstance(npts, (bool, np.bool_)) or npts < 3:
            raise ValueError("npts must be an integer >= 3")

        self.system = system
        self.potential = potential
        self.rmax = rmax
        self.npts = npts
        self.bc = bc

        # Uniform radial grid
        self.rgrid = np.linspace(0.0, rmax, npts)

        # C++ solver
        self._solver = rs.RadialSolver(
            system,
            potential,
            self.rgrid.tolist(),
            bc,
        )

    def solve(self, E):
        if not np.isfinite(E):
            raise ValueError("E must be finite")
        return self._solver.solve(E)