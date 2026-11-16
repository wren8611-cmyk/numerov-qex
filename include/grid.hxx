// grid.hxx

#pragma once

#include <vector>
#include <cstddef>
#include <stdexcept>

inline std::vector<double> linspace(double xmin, double xmax, std::size_t npts) {
    if (npts < 2)       throw std::invalid_argument("linspace requires at least two points");
    if (xmax <= xmin)   throw std::invalid_argument("xmax must be greater than xmin");

    std::vector<double> x(npts);
    
    const double h = (xmax - xmin) / (npts - 1);
    for (std::size_t i = 0; i < npts; ++i) x[i] = xmin + i * h;

    return x;
}