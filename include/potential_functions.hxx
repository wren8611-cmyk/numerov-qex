// potential_functions.hxx

#pragma once

#include <vector>
#include <string>
#include <stdexcept>
#include <cmath>

namespace potential_functions {

inline double evaluate(const std::string& name, const std::vector<double>& pars, double r) {
    if (name == "gaussian") {
        if (pars.size() != 2) throw std::invalid_argument("Gaussian requires 2 parameters");

        const double a = pars[0];
        const double b = pars[1];

        return a * std::exp(-std::pow(r / b, 2));
    }

    if (name == "well") {
        if (pars.size() != 2) throw std::invalid_argument("Well requires 2 parameters");

        const double a = pars[0];
        const double b = pars[1];

        return (r < b) ? a : 0.0;
    }

    throw std::invalid_argument("Unknown potential: " + name);
}

} // namespace potential_functions