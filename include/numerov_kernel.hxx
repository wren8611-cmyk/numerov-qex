// numerov_kernel.hxx

#pragma once

#include <Eigen/Dense>

#include <stdexcept>
#include <cstddef>

#include "channel_system.hxx"
#include "potential_matrix.hxx"

class NumerovKernel {
public:
    NumerovKernel(const ChannelSystem& system, const PotentialMatrix& potential) : system_(system), potential_(potential)
    { if (&potential_.channel_system() != &system_) throw std::invalid_argument("NumerovKernel: inconsistent ChannelSystem"); }

    Eigen::MatrixXcd eval_Q(double r, double E) const
    {
        constexpr double alpha = 1.0 / 137.035999084;
        constexpr double eps = 1e-14;

        const double r_reg = r + eps;

        const Eigen::Index n = static_cast<Eigen::Index>(system_.nchannels());

        Eigen::MatrixXcd Q = potential_.evaluate(r);

        // Q_ij = 2 mu_i V_ij
        for (Eigen::Index i = 0; i < n; ++i) {
            const auto& ch = system_[i];
            Q.row(i) *= 2.0 * ch.mu;
        }

        // Diagonal kinetic, centrifugal, and Coulomb terms
        for (Eigen::Index i = 0; i < n; ++i) {
            const auto&  ch  = system_[i];
            const double ell = static_cast<double>(ch.ell);

            Q(i, i) += ell * (ell + 1.0) / (r_reg * r_reg) + 2.0 * ch.mu * (ch.threshold - E + alpha * ch.charge_product / r_reg);
        }

        return Q;
    }

    const ChannelSystem& channel_system() const { return system_; }

protected:
    const ChannelSystem& system_;
    
private:
    const PotentialMatrix& potential_;
};