######################################################################
# sample_trajectory.jl
# computes a trajectory by sampling
######################################################################

using StatsBase: WeightVec, sample

# creates a sample_trajectory
function sample_trajectory(em::ErgodicManager, x0::T2F, h::Float64, N::Int)
	xd = Array(Vector{Float64}, N+1)
	xd[1] = [x0[1], x0[2]]
	ud = Array(Vector{Float64}, N+1)
	bin_size = (em.bins, em.bins)

	# first sample N points from e.phi
	weights = WeightVec(vec(em.phi))

	# if we don't use x0 and want slightly more randomness
	#xi, yi = ind2sub(bin_size, sample(weights))
	#xd[1] = [em.cell_size*(xi-.5), em.cell_size*(yi-.5)]

	for n = 1:N
		xi, yi = ind2sub(bin_size, sample(weights))
		xd[n+1] = [em.cell_size*(xi-.5), em.cell_size*(yi-.5)]
		ud[n] = (xd[n+1] - xd[n]) / h
	end
	ud[N+1] = ud[N]		# this one doesn't matter

	return xd, ud
end
