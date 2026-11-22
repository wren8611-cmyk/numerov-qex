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
    Calculate the partial-wave contribution Delta C_l(q).
    The full correlation function is
        C(q) = 1 + sum_l (2*l + 1) * Delta C_l(q)
    """
    if result.boundary_condition != BoundaryCondition.Incoming:
        raise ValueError("Incoming boundary condition required.")

    if not result.is_open(incoming_channel):
        raise ValueError("Incoming channel must be open.")

    channel = result.system.channel_by_id(incoming_channel)

    if channel.ell != ell:
        raise ValueError(f"Incoming channel has ell={channel.ell}, expected ell={ell}.")

    for channel_id in weights:
        ch = result.system.channel_by_id(channel_id)
        if ch.ell != ell:
            raise ValueError(f"Channel {channel_id} has ell={ch.ell}, expected ell={ell}.")

    q = np.sqrt(2.0 * channel.mu * (result.E - channel.threshold))

    if q <= 0:
        raise ValueError("Incoming momentum must be positive.")

    rgrid = np.asarray(result.wavefunction.x())
    r = rgrid[1:]

    norm = 2j * q
    correction = np.zeros_like(r, dtype=float)

    for channel_id, weight in weights.items():
        u = np.array([
            result.psi(i, channel_id, incoming_channel)
            for i in range(1, len(rgrid))
        ])

        correction += weight * np.abs(u / (norm * r))**2

    correction -= spherical_jn(ell, q * r)**2

    integrand = r**2 * source(r) * correction

    delta_C = 4.0 * np.pi * np.trapezoid(integrand, r)

    return q, float(delta_C)