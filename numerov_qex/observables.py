import numpy as np
from scipy.special import spherical_jn

from .radialsolver import BoundaryCondition


def phase_shift(result, channel_id):
    """
    Return the phase shift [rad] for a given channel.
    Returns nan if the channel is closed.
    """
    if result.boundary_condition != BoundaryCondition.Outgoing:
        raise ValueError("Outgoing boundary condition required.")

    sii = result.s_flux(channel_id, channel_id)

    return 0.5 * np.angle(sii)


def inelasticity(result, channel_id):
    """
    Return the inelasticity for a given channel.
    Returns nan if the channel is closed.
    """
    if result.boundary_condition != BoundaryCondition.Outgoing:
        raise ValueError("Outgoing boundary condition required.")

    sii = result.s_flux(channel_id, channel_id)

    return np.abs(sii)


def correlation_function(result, source, weights, incoming_channel, ell):
    """
    Calculate the partial-wave femtoscopic correlation function C_l(q).

    Parameters
    ----------
    result : ScatteringResult
        Scattering solution.
    source : callable
        Normalized radial source function S(r).
    weights : dict
        Channel weights {channel_id: weight}.
    incoming_channel : int
        Incoming channel ID.
    ell : int
        Orbital angular momentum.

    Returns
    -------
    q : float
        Relative momentum in the incoming channel.
    C : float
        Partial-wave correlation function C_l(q).
    """
    if result.boundary_condition != BoundaryCondition.Incoming:
        raise ValueError("Incoming boundary condition required.")

    if not result.is_open(incoming_channel):
        raise ValueError("Incoming channel must be open.")

    # Check orbital angular momenta
    channel = result.system.channel_by_id(incoming_channel)

    if channel.ell != ell:
        raise ValueError(f"Incoming channel has ell={channel.ell}, expected ell={ell}.")

    for channel_id in weights:
        ch = result.system.channel_by_id(channel_id)
        if ch.ell != ell:
            raise ValueError(f"Channel {channel_id} has ell={ch.ell}, expected ell={ell}.")

    # Incoming momentum
    q = np.sqrt(2.0 * channel.mu * (result.E - channel.threshold))

    if q <= 0:
        raise ValueError("Incoming momentum must be positive.")

    # Radial grid
    rgrid = np.asarray(result.wavefunction.x())
    r = rgrid[1:]

    norm = 2j * q
    correction = np.zeros_like(r, dtype=float)

    # Coupled-channel wave functions
    for channel_id, weight in weights.items():
        u = np.array([result.psi(i, channel_id, incoming_channel) for i in range(1, len(rgrid))])

        correction += weight * np.abs(u / (norm * r))**2

    # Subtract free partial wave
    correction -= spherical_jn(ell, q * r)**2

    # Source integration
    integrand = r**2 * source(r) * correction

    return q, float(1.0 + 4.0 * np.pi * np.trapezoid(integrand, r))