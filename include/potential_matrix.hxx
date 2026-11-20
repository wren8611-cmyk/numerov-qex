// potential_matrix.hxx

#pragma once

#include "channel_system.hxx"
#include "potential_info.hxx"
#include "potential_functions.hxx"

#include <Eigen/Dense>

#include <iostream>
#include <vector>

class PotentialMatrix {
public:
    using Matrix = Eigen::MatrixXcd;

    explicit PotentialMatrix(const std::vector<PotentialInfo>& potentials, const ChannelSystem& system)
    : system_(system) { append(potentials); }

    void append(const PotentialInfo& potential) {
        // Validate channel IDs.
        system_.channel_index(potential.channelid_left);
        system_.channel_index(potential.channelid_right);
        potentials_.push_back(potential);
    }
    void append(const std::vector<PotentialInfo>& potentials) { for (const auto& p : potentials) append(p); }

    void clear() { potentials_.clear(); }

    void display_info(std::ostream& os = std::cout) const {
        if (empty()){
            os << "POTENTIAL MATRIX EMPTY" << "\n";
            return;
        }
        
        for (const auto& p : potentials_) p.summary(os);
    }

    std::size_t size() const {
        return potentials_.size();
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

    bool empty() const {
        return potentials_.empty();
    }

};