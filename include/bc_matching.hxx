// bc_matching.hxx
#pragma once

#include "channel_system.hxx"
#include "wavefunction.hxx"
#include "coulomb_functions.hxx"
#include "closed_coulomb_functions.hxx"

#include <Eigen/Dense>

#include <algorithm>
#include <vector>
#include <complex>
#include <cmath>
#include <stdexcept>
#include <cstddef>
#include <limits>
#include <utility>

enum class BoundaryCondition {
    Outgoing,
    Incoming
};

struct ScatteringResult {
    using Matrix = Eigen::MatrixXcd;
    using cd = std::complex<double>;

    Matrix S;
    Matrix S_flux;
    MatrixFunction wavefunction;

    double E;
    const ChannelSystem& system;

    std::vector<int> open_channel_ids;
    BoundaryCondition boundary_condition;

    // Convert open-channel ID to matrix column index.
    std::size_t open_channel_index(int channel_id) const {
        for (std::size_t a = 0; a < open_channel_ids.size(); ++a) {
            if (open_channel_ids[a] == channel_id) return a;
        }

        throw std::out_of_range("Channel ID is not an open channel");
    }

    // Physical scattering wavefunction component.
    cd psi(std::size_t radial_index, int channel_id, int reference_channel_id) const {
        return wavefunction[radial_index](
            system.channel_index(channel_id),
            open_channel_index(reference_channel_id)
        );
    }

    bool is_open(int channel_id) const {
        system.channel_index(channel_id); // Validate channel ID.
        return std::find( open_channel_ids.begin(), open_channel_ids.end(), channel_id) != open_channel_ids.end();
    }

    // accessor
    cd s(int channel_id_left, int channel_id_right) const {
        const double nan = std::numeric_limits<double>::quiet_NaN();

        if (!is_open(channel_id_left) || !is_open(channel_id_right)) return cd(nan, nan);
        return S(open_channel_index(channel_id_left), open_channel_index(channel_id_right));
    }

    cd s_flux(int channel_id_left, int channel_id_right) const {
        const double nan = std::numeric_limits<double>::quiet_NaN();

        if (!is_open(channel_id_left) || !is_open(channel_id_right)) return cd(nan, nan);
        return S_flux(open_channel_index(channel_id_left), open_channel_index(channel_id_right));
    }

};

class BCMatching {
public:
    using Matrix = Eigen::MatrixXcd;
    using cd = std::complex<double>;

    explicit BCMatching(const ChannelSystem& system, BoundaryCondition bc = BoundaryCondition::Outgoing)
        : system_(system), bc_(bc) {}

    ScatteringResult match(const Wavefunction& wf) const
    {

        if (&wf.channel_system() != &system_)
            throw std::invalid_argument("BCMatching: inconsistent ChannelSystem");
            
        const double E = wf.energy();
        const std::size_t N = system_.nchannels();
        const std::size_t M = wf.size();

        if (M < 3)
            throw std::invalid_argument(
                "At least three radial points are required"
            );

        const auto& radius = wf.x();

        const std::size_t n = M - 2;

        const double R = radius[n];
        const double h = radius[1] - radius[0];

        if (!std::isfinite(R) || !std::isfinite(h) ||
            R <= 0.0 || h <= 0.0)
            throw std::invalid_argument(
                "Invalid matching radius or grid spacing"
            );

        constexpr double tol_abs = 1e-12;
        constexpr double tol_rel = 1e-8;

        for (std::size_t i = 1; i < M; ++i) {
            const double hi = radius[i] - radius[i - 1];

            if (!std::isfinite(hi) || hi <= 0.0 ||
                std::abs(hi - h) > std::max(tol_abs, tol_rel * h))
                throw std::invalid_argument(
                    "Radial grid must be uniform and strictly increasing"
                );
        }

        // Numerov wavefunction and its radial derivative.
        const Matrix U = wf[n];

        const Matrix Uprime =
            (wf[n + 1] - wf[n - 1]) / (2.0 * h);

        // Identify open and closed channels.
        std::vector<std::size_t> open;
        std::vector<std::size_t> closed;

        for (std::size_t i = 0; i < N; ++i) {
            const auto& ch = system_[i];

            if (E > ch.threshold) {
                open.push_back(i);
            }
            else if (E < ch.threshold) {
                closed.push_back(i);
            }
            else {
                throw std::invalid_argument(
                    "Matching exactly at threshold is not supported"
                );
            }
        }

        const std::size_t No = open.size();

        if (No == 0)
            throw std::invalid_argument(
                "At least one open channel is required"
            );

        // Incoming asymptotic wavefunctions.
        Matrix I = Matrix::Zero(N, No);
        Matrix Iprime = Matrix::Zero(N, No);

        // Outgoing/decaying reference functions
        // and their logarithmic derivatives.
        Eigen::VectorXcd O(N);
        Eigen::VectorXcd L(N);

        // Open channels: Coulomb functions.
        for (std::size_t a = 0; a < No; ++a) {
            const std::size_t i = open[a];
            const auto& ch = system_[i];

            const double k =
                std::sqrt(2.0 * ch.mu * (E - ch.threshold));

            const double eta =
                alpha * ch.mu * ch.charge_product / k;

            const auto cf =
                coulomb_functions(ch.ell, eta, k, R);

            if (bc_ == BoundaryCondition::Outgoing) {
                O(i) = cf.Hplus;
                L(i) = cf.dHplus / cf.Hplus;

                I(i, a) = cf.Hminus;
                Iprime(i, a) = cf.dHminus;
            }
            else if (bc_ == BoundaryCondition::Incoming) {
                O(i) = cf.Hminus;
                L(i) = cf.dHminus / cf.Hminus;

                I(i, a) = cf.Hplus;
                Iprime(i, a) = cf.dHplus;
            }
            else { throw std::invalid_argument("Invalid boundary condition"); }
        }

        // Closed channels: decaying Coulomb wavefunctions.
        for (const auto i : closed) {
            const auto& ch = system_[i];

            const double kappa =
                std::sqrt(2.0 * ch.mu * (ch.threshold - E));

            const double eta =
                alpha * ch.mu * ch.charge_product / kappa;

            O(i) = 1.0;

            L(i) = closed_coulomb_logderivative(
                ch.ell, eta, kappa, R
            );
        }

        // Solve:
        //
        // (U' - L U) C = I' - L I
        //
        // L acts on the channel (row) index.

        Matrix A = Uprime;
        Matrix B = Iprime;

        for (std::size_t i = 0; i < N; ++i) {
            A.row(i) -= L(i) * U.row(i);
            B.row(i) -= L(i) * I.row(i);
        }

        const Matrix C = A.partialPivLu().solve(B);
        
        const Matrix UC = U * C;

        Matrix S(No, No);

        for (std::size_t a = 0; a < No; ++a) {
            const std::size_t i = open[a];
            S.row(a) = (I.row(i) - UC.row(i)) / O(i);
        }

        // Flux-normalized scattering matrix.
        Eigen::VectorXd sqrt_v(No);

        for (std::size_t a = 0; a < No; ++a) {
            const auto& ch = system_[open[a]];

            const double k =
                std::sqrt(2.0 * ch.mu * (E - ch.threshold));

            sqrt_v(a) = std::sqrt(k / ch.mu);
        }

        Matrix S_flux = S;

        for (std::size_t a = 0; a < No; ++a) {
            for (std::size_t b = 0; b < No; ++b) {
                S_flux(a, b) *= sqrt_v(a) / sqrt_v(b);
            }
        }

        // Transform regular solutions into physical scattering solutions.
        MatrixFunction physical_wf = wf.matmul_right(C);

        std::vector<int> open_channel_ids;
        open_channel_ids.reserve(No);
        for (const auto i : open) open_channel_ids.push_back(system_[i].id);

        return {
            std::move(S),
            std::move(S_flux),
            std::move(physical_wf),
            E,
            system_,
            std::move(open_channel_ids),
            bc_
        };
    }

private:
    const ChannelSystem& system_;
    BoundaryCondition bc_;

    static constexpr double alpha = 1.0 / 137.035999084;
};