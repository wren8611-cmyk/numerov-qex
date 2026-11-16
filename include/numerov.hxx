#pragma once

#include "channel_system.hxx"
#include "potential_matrix.hxx"
#include "wavefunction.hxx"
#include "numerov_kernel.hxx"

#include <Eigen/Dense>

#include <algorithm>
#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <stdexcept>
#include <cstddef>
#include <utility>

class Numerov : public NumerovKernel {
public:
    using Matrix = Eigen::MatrixXcd;

    Numerov(const ChannelSystem& system, const PotentialMatrix& potential, std::vector<double> radius)
        : NumerovKernel(system, potential), radius_(std::move(radius))
    {
        validate_grid();
    }

    // Integrate regular solutions.
    Wavefunction solve(double E) const {
        const std::size_t N = system_.nchannels();
        const std::size_t M = radius_.size();

        const double h = radius_[1] - radius_[0];

        Wavefunction wf(radius_, system_, E);

        // Regular boundary conditions.
        wf[0] = Matrix::Zero(N, N);
        wf[1] = Matrix::Zero(N, N);

        for (std::size_t i = 0; i < N; ++i) {
            const int ell = system_[i].ell;
            wf[1](i, i) = std::pow(h, ell + 1);
        }

        const Matrix I = Matrix::Identity(N, N);
        const double h2 = h * h;

        Matrix Q_prev = eval_Q(radius_[0], E);
        Matrix Q_curr = eval_Q(radius_[1], E);

        // Matrix Numerov recurrence.
        constexpr double warning_threshold = 1e100;
        bool warned = false;
        for (std::size_t n = 1; n + 1 < M; ++n) {
            const Matrix Q_next = eval_Q(radius_[n + 1], E);

            const Matrix A_prev = I - (h2 / 12.0) * Q_prev;
            const Matrix B_curr = I + (5.0 * h2 / 12.0) * Q_curr;
            const Matrix A_next = I - (h2 / 12.0) * Q_next;

            const Matrix rhs = 2.0 * B_curr * wf[n] - A_prev * wf[n - 1];

            wf[n + 1] = A_next.partialPivLu().solve(rhs);

            const double max_abs = wf[n + 1].cwiseAbs().maxCoeff();

            if (!warned && (!std::isfinite(max_abs) || max_abs > warning_threshold)) {
                std::cerr
                    << "Warning: Numerov wavefunction is excessively large"
                    << " at r = " << radius_[n + 1]
                    << ", max|U| = " << max_abs << '\n';

                warned = true;
            }

            Q_prev = Q_curr;
            Q_curr = Q_next;
        }

        return wf;
    }

private:
    std::vector<double> radius_;

    void validate_grid() const {
        const std::size_t M = radius_.size();

        if (M < 3) throw std::invalid_argument("At least three radial points are required");

        if (!std::isfinite(radius_[0]) || std::abs(radius_[0]) > 1e-15)
            throw std::invalid_argument("Radial grid must start at r = 0");

        const double h = radius_[1] - radius_[0];

        if (!std::isfinite(h) || h <= 0.0) throw std::invalid_argument("Radial grid spacing must be positive");

        constexpr double tol_abs = 1e-12;
        constexpr double tol_rel = 1e-8;
        for (std::size_t i = 1; i < M; ++i) {
            const double hi = radius_[i] - radius_[i - 1];

            if (!std::isfinite(hi) || hi <= 0.0 || std::abs(hi - h) > std::max(tol_abs, tol_rel * h))
                throw std::invalid_argument("Numerov requires a uniform radial grid");
        }
    }
};