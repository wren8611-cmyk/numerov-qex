// coulomb_functions.hxx

#pragma once

#include <gsl/gsl_sf_coulomb.h>
#include <gsl/gsl_errno.h>

#include <complex>
#include <cmath>
#include <stdexcept>

struct CoulombFunctions {
    std::complex<double> Hplus;
    std::complex<double> Hminus;
    std::complex<double> dHplus;
    std::complex<double> dHminus;
};

inline CoulombFunctions coulomb_functions(
    int ell,
    double eta,
    double k,
    double r
) {
    if (ell < 0 || k <= 0.0 || r <= 0.0)
        throw std::invalid_argument(
            "Invalid Coulomb function arguments"
        );

    const double rho = k * r;

    gsl_sf_result F, Fp, G, Gp;
    double exp_F, exp_G;

    const int status = gsl_sf_coulomb_wave_FG_e(
        eta,
        rho,
        static_cast<double>(ell),
        0,
        &F,
        &Fp,
        &G,
        &Gp,
        &exp_F,
        &exp_G
    );

    if (status != GSL_SUCCESS)
        throw std::runtime_error(
            "GSL Coulomb function evaluation failed"
        );

    // GSL may return scaled values when the functions
    // are too large for ordinary double precision.
    const double scale_F = std::exp(exp_F);
    const double scale_G = std::exp(exp_G);

    const double f  = F.val  * scale_F;
    const double fp = Fp.val * scale_F;
    const double g  = G.val  * scale_G;
    const double gp = Gp.val * scale_G;

    const std::complex<double> imag(0.0, 1.0);

    return {
        g + imag * f,
        g - imag * f,
        k * (gp + imag * fp),
        k * (gp - imag * fp)
    };
}