// wavefunction.hxx

#pragma once

#include "matrix_function.hxx"
#include "channel_system.hxx"

class Wavefunction : public MatrixFunction {
public:
    // Constructor
    Wavefunction(const std::vector<double>& radius, const ChannelSystem& system, const double E)
        : MatrixFunction(radius, system.nchannels()), system_(system), E_(E) {}

    // Accessors
    const ChannelSystem& channel_system() const { return system_; }
    double energy() const { return E_; }

private:
    const ChannelSystem& system_;
    const double E_;
};