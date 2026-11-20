#pragma once

#include <vector>
#include <cstddef>
#include <stdexcept>
#include <cmath>
#include <algorithm>
#include <utility>

struct RadialGrid {
public:
    explicit RadialGrid(std::vector<double> radius) : radius_(std::move(radius)) {
        validate();
        h_ = radius_[1] - radius_[0];
    }

    RadialGrid(double rmax, std::size_t npts) : radius_(npts) {
        if (!std::isfinite(rmax) || rmax <= 0.0)    throw std::invalid_argument("rmax must be positive and finite");
        if (npts < 3)                               throw std::invalid_argument("At least three radial points are required");

        h_ = rmax / (npts - 1);
        for (std::size_t i = 0; i < npts; ++i) radius_[i] = i * h_;

        validate();
    }

    // Accessors
    const std::vector<double>& values() const { return radius_; }
    double operator[](std::size_t i)    const { return radius_[i]; }
    double h()                          const { return h_; }
    std::size_t size()                  const { return radius_.size(); }
    double rmax()                       const { return radius_.back(); }

private:
    std::vector<double> radius_;
    double h_;

    void validate() const {
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