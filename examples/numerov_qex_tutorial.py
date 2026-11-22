"""
Tutorial for numerov-qex solver

Use consistent natural units (i.e. hbar = c = 1). 
Mass, energy, momentum, and potential strength share one unit. 
"""
import numpy as np
import numerov_qex.radialsolver as rs
import numerov_qex.solver as nv
from numerov_qex.observables import inelasticity, phase_shift


# Specify physical system
# ChannelInfo(id, name, mass_1, mass_2, orbital_angular_momentum, charge_1, charge_2)
# id and name can be anything convenient
CHANNEL0_ID = 4649
CHANNEL1_ID = 67
system = rs.ChannelSystem([
    rs.ChannelInfo(CHANNEL0_ID, "channel0 in s-wave", 1.0, 1.0, 0, 0.0, 0.0),
    rs.ChannelInfo(CHANNEL1_ID, "channel1 in s-wave", 1.0, 1.5, 0, 0.0, 0.0),
])
system.info()
print("\n\n\n\n\n\n")



# Set up the potential matrix
# Leave it empty for free motion in neutral channels. 
# If the charge product is nonzero, the solver still includes diagonal Coulomb (automatic).
potential = rs.PotentialMatrix([], system)

# Append potentials
# append() accumulates terms, even at the same matrix element.
# To start a new potential, call potential.clear() before appending its terms.
# The solver sees changes to this PotentialMatrix at evaluation time not solver-construction time.
# Coulomb is auto-set.

# currently support: "gaussian", "well"
# definitions in `(base)/include/potential_functions.hxx`
# PotentialInfo("gaussian", [a, b], row_channel_id, column_channel_id). ### a * exp[-(r / b)^2]
potential.append(rs.PotentialInfo("gaussian", [-5.0, 1.0], CHANNEL0_ID, CHANNEL0_ID))
potential.append([
    rs.PotentialInfo("gaussian", [-3.0, 1.0], CHANNEL1_ID, CHANNEL1_ID),
    rs.PotentialInfo("gaussian", [-1.0, 1.0], CHANNEL0_ID, CHANNEL1_ID),
    rs.PotentialInfo("gaussian", [-1.0, 1.0], CHANNEL1_ID, CHANNEL0_ID),
])  # Add both off-diagonal entries (the solver will not check hermiticity)
potential.display_info()
print(f"Stored terms: {potential.size()}")
print("Explicit potential matrix at r = 1:\n", potential.evaluate(1.0))
print("\n\n\n\n\n\n")




# Set up solver
# supported bc: "BoundaryCondition.Outgoing", "BoundaryCondition.Incoming"
solver = nv.NumerovSolver(
    system, potential, 
    rmax=20.0, 
    npts=2048,
    bc=rs.BoundaryCondition.Outgoing,
)





# Solve at CM-energy E (including the rest masses)
# Both channels are open at 2.7. Try 2.2 for closed channel 1.
result = solver.solve(2.7)
#result = solver.solve(2.2)





# Result object
print(f"E = {result.E:g}, boundary condition = {result.boundary_condition}")
print("Open channel IDs in S-matrix order:", result.open_channel_ids)

for channel_id in (CHANNEL0_ID, CHANNEL1_ID):
    channel = result.system.channel_by_id(channel_id)
    print(
        f"{channel.name}: mu = {channel.mu:g}, threshold = {channel.threshold:g}, "
        f"ell = {channel.ell}, charge product = {channel.charge_product:g}"
    )
    if result.is_open(channel_id):
        index = result.open_channel_index(channel_id)
        delta = phase_shift(result, channel_id)
        eta   = inelasticity(result, channel_id)
        print(f"  S-matrix index {index}, delta = {delta:+.6f} rad, eta = {eta:.6f}")
    else:
        print("  Closed")

print("S matrix:\n", result.S)
print("Flux-normalized S matrix:\n", result.S_flux)
print("det(S^dagger S): ", np.linalg.det((result.S).conj().T @ (result.S)))
print("S_flux^dagger S_flux:\n", (result.S_flux).conj().T @ (result.S_flux))
# The scalar accessors result.s() take channel IDs, while arrays use open-channel indices (result.S[i,j] != result.s(i,j))
print("S(0, 0):", result.s(CHANNEL0_ID, CHANNEL0_ID))
print("S(1, 1):", result.s(CHANNEL1_ID, CHANNEL1_ID))
print("S_flux(0, 0):", result.s_flux(CHANNEL0_ID, CHANNEL0_ID))
print("S_flux(1, 1):", result.s_flux(CHANNEL1_ID, CHANNEL1_ID))
print("\n")

wavefunction = result.wavefunction
radii = wavefunction.x()
radial_index = len(wavefunction) // 2
print(f"Stored wavefunction: {len(wavefunction)} radial points")
print(f"Matrix shape at each point: {wavefunction.nrows()} x {wavefunction.ncols()}")
print(f"Wavefunction matrix at r = {radii[radial_index]:.4f}:\n", wavefunction[radial_index])
# psi(radial_index, component_channel_id, reference_channel_id).
# The reference must be open. The component can be open or closed.
# Note that its specified by channel_id
print("Channel 1 component for reference channel 0:", result.psi(radial_index, CHANNEL1_ID, CHANNEL0_ID))
print("\n\n\n\n\n\n")

# To calculate a fresh potential with the same solver later:
# potential.clear()
# potential.append(rs.PotentialInfo("gaussian", [-2.0, 1.0], CHANNEL0_ID, CHANNEL0_ID))
# updated_result = solver.solve(result.E)