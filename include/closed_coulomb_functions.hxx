// closed_coulomb_functions.hxx

#pragma once

#include <cmath>
#include <stdexcept>
#include <algorithm>

inline double closed_coulomb_logderivative(
    int ell,
    double eta,
    double kappa,
    double r,
    double rmax = 0.0,
    double dr = 0.01
) {
    if (ell < 0)
        throw std::invalid_argument("ell must be nonnegative");

    if (kappa <= 0.0 || r <= 0.0 || dr <= 0.0)
        throw std::invalid_argument("kappa, r, and dr must be positive");

    // Choose a sufficiently large starting radius.
    if (rmax <= 0.0)
        rmax = r + std::max(20.0 / kappa, 20.0);

    if (rmax <= r)
        throw std::invalid_argument("rmax must be greater than r");

    // Riccati equation:
    //
    // L' = Q(r) - L^2
    //
    // Q(r) = kappa^2 + 2*kappa*eta/r + ell*(ell+1)/r^2

    auto rhs = [&](double x, double L) {
        const double Q =
            kappa * kappa
            + 2.0 * kappa * eta / x
            + static_cast<double>(ell * (ell + 1)) / (x * x);

        return Q - L * L;
    };

    // Asymptotic boundary condition.
    //
    // W_{-eta,ell+1/2}(2*kappa*r)
    // ~ exp(-kappa*r) * (2*kappa*r)^(-eta)
    //
    // Therefore L ~ -kappa - eta/r.

    double L = -kappa - eta / rmax;

    const int nsteps =
        std::max(1, static_cast<int>(std::ceil((rmax - r) / dr)));

    const double h = (r - rmax) / nsteps;

    double x = rmax;

    for (int n = 0; n < nsteps; ++n) {
        const double k1 = rhs(x, L);
        const double k2 = rhs(x + 0.5 * h, L + 0.5 * h * k1);
        const double k3 = rhs(x + 0.5 * h, L + 0.5 * h * k2);
        const double k4 = rhs(x + h, L + h * k3);

        L += (h / 6.0) * (k1 + 2.0*k2 + 2.0*k3 + k4);
        x += h;
    }

    return L;
}