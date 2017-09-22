######################################################################
# ergodic_manager.jl
######################################################################
abstract type ErgodicManager end

decompose!(em::ErgodicManager) = decompose!(em, em.phi)

"""
`phi = reconstruct(em)`

Reconstructs distribution from Fourier coefficients `em.phik`.

`phi = reconstruct(em, ck::Matrix{Float64})`

Reconstructs distribution from Fourier coefficients `ck`.

Both options return a distribution `phi`, a `em.bins` by `em.bins` matrix.
The value at `phi[xi,yi]` gives the value at the cell index `xi,yi`.

`phi = reconstruct(em, xd::VVF)`

Reconstructs from trajectory `xd`.
"""
reconstruct(em::ErgodicManager) = reconstruct(em, em.phik)
function reconstruct(em::ErgodicManager, xd::VVF)
	reconstruct(em, decompose(em, xd))
end


x_min(em::ErgodicManager) = em.domain.mins[1]
y_min(em::ErgodicManager) = em.domain.mins[2]
z_min(em::ErgodicManager) = em.domain.mins[3]

x_max(em::ErgodicManager) = em.domain.maxes[1]
y_max(em::ErgodicManager) = em.domain.maxes[2]
z_max(em::ErgodicManager) = em.domain.maxes[3]

x_size(em::ErgodicManager) = em.domain.cell_lengths[1]
y_size(em::ErgodicManager) = em.domain.cell_lengths[2]
z_size(em::ErgodicManager) = em.domain.cell_lengths[3]

x_cells(em::ErgodicManager) = em.domain.cells[1]
y_cells(em::ErgodicManager) = em.domain.cells[2]
z_cells(em::ErgodicManager) = em.domain.cells[3]


# decomposition on a group of trajectories
function decompose(em::ErgodicManager, xds::VVVF)
	num_agents = length(xds)

	ck = decompose(em, xds[1])
	for i = 2:num_agents
		ck += decompose(em, xds[i])
	end

	# divide by the number of agents (that is, length(xds))
	return ck / num_agents
end


"""
`ergodic_score(em, traj::VVF)`

First breaks down the trajectory into components ck.
"""
function ergodic_score(em::ErgodicManager, traj::VVF)
	ck = decompose(em, traj)
	return ergodic_score(em, ck)
end


function ergodic_score(em::ErgodicManager, ck)
	val = 0.0
	for (i, L) in enumerate(em.Lambda)
		d = em.phik[i] - ck[i]
		val += L * d * d
	end
	return val
end

normalize!(em::ErgodicManager) = normalize!(em.phi, em.domain.cell_size)
