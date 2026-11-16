// potential_matrix.hxx

#pragma once

#include "channel_system.hxx"
#include "potential_info.hxx"
#include "potential_functions.hxx"

#include <Eigen/Dense>

#include <vector>
#include <stdexcept>

class PotentialMatrix {
public:
    using Matrix = Eigen::MatrixXcd;

    PotentialMatrix(const std::vector<PotentialInfo>& potentials, const ChannelSystem& system)
        : potentials_(potentials), system_(system)
    {
        for (const auto& p : potentials_) {
            system_.channel_index(p.channelid_left);
            system_.channel_index(p.channelid_right);
        }
    }

    Matrix evaluate(double r) const {
        Matrix V = Matrix::Zero(
            system_.nchannels(),
            system_.nchannels()
        );

        for (const auto& p : potentials_) {
            const auto i = system_.channel_index(p.channelid_left);
            const auto j = system_.channel_index(p.channelid_right);

            V(i, j) += potential_functions::evaluate(p.name, p.pars, r);
        }

        return V;
    }

    const ChannelSystem& channel_system() const { return system_; }

private:
    std::vector<PotentialInfo> potentials_;
    const ChannelSystem& system_;
};