// radial_solver.hxx

#pragma once

#include "channel_system.hxx"
#include "potential_matrix.hxx"
#include "numerov.hxx"
#include "bc_matching.hxx"

#include <vector>

class RadialSolver {
public:
    // Constructors
    RadialSolver(
        const ChannelSystem& system,
        const PotentialMatrix& potential,
        const RadialGrid& rgrid,
        BoundaryCondition bc = BoundaryCondition::Outgoing
    ) : numerov_(system, potential, rgrid), matching_(system, bc) {}
    RadialSolver(
        const ChannelSystem& system,
        const PotentialMatrix& potential,
        const std::vector<double>& radius,
        BoundaryCondition bc = BoundaryCondition::Outgoing
    ) : RadialSolver(system, potential, RadialGrid(radius), bc) {}

    // Accessor
    const RadialGrid& rgrid() const { return numerov_.rgrid(); }

    // Solve the coupled-channel Schrodinger equation and match asymptotics.
    ScatteringResult solve(double E) const {
        auto wf = numerov_.solve(E);
        return matching_.match(wf);
    }

private:
    Numerov numerov_;
    BCMatching matching_;
};