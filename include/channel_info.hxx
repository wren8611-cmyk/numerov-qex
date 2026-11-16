// channel_info.hxx 

#pragma once

#include <iomanip>
#include <iostream>
#include <string>
#include <stdexcept>
#include <utility>

struct ChannelInfo {
    int id;
    std::string name;
    double mu;              // reduced mass
    double threshold;       // threshold mass
    int ell;                // orb angular momentum
    double charge_product;

    ChannelInfo(int id, std::string name, double m1, double m2, int ell = 0, double z1 = 0.0, double z2 = 0.0)
        : id(id), name(std::move(name)), mu(m1 * m2 / (m1 + m2)), threshold(m1 + m2), ell(ell), charge_product(z1 * z2)
    {
        if (m1 <= 0.0 || m2 <= 0.0) throw std::invalid_argument("Particle masses must be positive");
        if (ell < 0)                throw std::invalid_argument("Orbital angular momentum must be nonnegative");
    }

    // Display channel information
    void summary(std::ostream& os = std::cout) const {
        os << "Channel Information\n";
        os << "-------------------\n";
        os << "ID             : " << id << '\n';
        os << "Name           : " << name << '\n';
        os << "Reduced mass   : " << mu << '\n';
        os << "Threshold      : " << threshold << '\n';
        os << "ell            : " << ell << '\n';
        os << "Charge product : " << charge_product << '\n';
    }
};