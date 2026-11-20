// wavefunction.hxx

#pragma once

#include "matrix_function.hxx"
#include "channel_system.hxx"

class Wavefunction : public MatrixFunction {
public:
    // Constructors
    Wavefunction(const RadialGrid& rgrid, const ChannelSystem& system, double E)
        : MatrixFunction(rgrid, system.nchannels()), system_(system), E_(E) {}

    Wavefunction(const std::vector<double>& radius, const ChannelSystem& system, double E)
        : Wavefunction(RadialGrid(radius), system, E) {}

    // Accessors
    const ChannelSystem& channel_system() const { return system_; }
    double energy() const { return E_; }

private:
    const ChannelSystem& system_;
    const double E_;
};