// radial_solver.hxx

#pragma once

#include "channel_system.hxx"
#include "potential_matrix.hxx"
#include "numerov.hxx"
#include "bc_matching.hxx"

#include <vector>
#include <utility>

class RadialSolver {
public:
    RadialSolver(const ChannelSystem& system, const PotentialMatrix& potential, std::vector<double> rgrid, const BoundaryCondition bc = BoundaryCondition::Outgoing)
        : numerov_(system, potential, std::move(rgrid)), matching_(system, bc)
    {}

    // Solve the coupled-channel Schrodinger equation and match asymptotics.
    ScatteringResult solve(const double E) const {
        auto wf = numerov_.solve(E);
        return matching_.match(wf);
    }

private:
    Numerov numerov_;
    BCMatching matching_;
};