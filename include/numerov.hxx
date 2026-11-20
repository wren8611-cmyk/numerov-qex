#pragma once

#include "channel_system.hxx"
#include "potential_matrix.hxx"
#include "wavefunction.hxx"
#include "numerov_kernel.hxx"
#include "grid.hxx"

#include <Eigen/Dense>

#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <cstddef>

class Numerov : public NumerovKernel {
public:
    using Matrix = Eigen::MatrixXcd;

    // Constructors
    Numerov(const ChannelSystem& system, const PotentialMatrix& potential, const RadialGrid& rgrid)
        : NumerovKernel(system, potential), rgrid_(rgrid) {}
    Numerov(const ChannelSystem& system, const PotentialMatrix& potential, const std::vector<double>& radius)
        : Numerov(system, potential, RadialGrid(radius)) {}

    // Accessor
    const RadialGrid& rgrid() const { return rgrid_; }

    // Integrate regular solutions.
    Wavefunction solve(double E) const {
        const std::size_t N = system_.nchannels();
        const std::size_t M = rgrid_.size();

        const double h = rgrid_.h();

        Wavefunction wf(rgrid_, system_, E);

        // Regular boundary conditions.
        wf[0] = Matrix::Zero(N, N);
        wf[1] = Matrix::Zero(N, N);

        for (std::size_t i = 0; i < N; ++i) {
            const int ell = system_[i].ell;
            wf[1](i, i) = std::pow(h, ell + 1);
        }

        const Matrix I = Matrix::Identity(N, N);
        const double h2 = h * h;

        Matrix Q_prev = eval_Q(rgrid_[0], E);
        Matrix Q_curr = eval_Q(rgrid_[1], E);

        // Matrix Numerov recurrence.
        constexpr double warning_threshold = 1e100;
        bool warned = false;

        for (std::size_t n = 1; n + 1 < M; ++n) {
            const Matrix Q_next = eval_Q(rgrid_[n + 1], E);

            const Matrix A_prev = I - (h2 / 12.0) * Q_prev;
            const Matrix B_curr = I + (5.0 * h2 / 12.0) * Q_curr;
            const Matrix A_next = I - (h2 / 12.0) * Q_next;

            const Matrix rhs = 2.0 * B_curr * wf[n] - A_prev * wf[n - 1];

            wf[n + 1] = A_next.partialPivLu().solve(rhs);

            const double max_abs = wf[n + 1].cwiseAbs().maxCoeff();

            if (!warned && (!std::isfinite(max_abs) || max_abs > warning_threshold)) {
                std::cerr
                    << "Warning: Numerov wavefunction is excessively large"
                    << " at r = " << rgrid_[n + 1]
                    << ", max|U| = " << max_abs << '\n';

                warned = true;
            }

            Q_prev = Q_curr;
            Q_curr = Q_next;
        }

        return wf;
    }

private:
    RadialGrid rgrid_;
};