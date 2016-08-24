######################################################################
# sample_trajectory.jl
# computes a trajectory by sampling
######################################################################

using StatsBase: WeightVec, sample

# creates a sample_trajectory
function sample_trajectory(em::ErgodicManager, tm::TrajectoryManager)
	x0 = (tm.x0[1], tm.x0[2])
	return sample_trajectory(em, x0, tm.h, tm.N)
end
function sample_trajectory(em::ErgodicManager, x0::T2F, h::Float64, N::Int)
	xd = Array(Vector{Float64}, N+1)
	points = Array(Vector{Float64}, N)
	xd[1] = [x0[1], x0[2]]
	ud = Array(Vector{Float64}, N)
	bin_size = (em.bins, em.bins)

	# first sample N points from e.phi
	weights = WeightVec(vec(em.phi))
	for n = 1:N
		xi, yi = ind2sub(bin_size, sample(weights))
		points[n] = [em.cell_size*(xi-.5), em.cell_size*(yi-.5)]
	end

	# find a short path heuristically
	#tsp_rand!(xd, points)
	tsp_nn!(xd, points)

	# compute controls
	ud = compute_controls(xd, h)

	return xd, ud
end



# simply takes the points as they are
# this creates a shitty trajectory
function tsp_rand!(xd::VV_F, points::VV_F)
	n = 1
	for p in points
		xd[n+1] = p
		n += 1
	end
end
