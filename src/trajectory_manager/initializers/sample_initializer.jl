######################################################################
# sample_initializer.jl
######################################################################

"""
`si = SampleInitializer()`

Samples points from a distribution.
"""
mutable struct SampleInitializer <: Initializer end


function initialize(si::SampleInitializer, em::ErgodicManager, tm::TrajectoryManager)
	xd = Array(Vector{Float64}, tm.N+1)
	points = Array(Vector{Float64}, tm.N)
	xd[1] = deepcopy(tm.x0)
	ud = Array(Vector{Float64}, tm.N)
	bin_size = (em.bins, em.bins)

	# first sample N points from e.phi
	weights = Weights(vec(em.phi))
	for n = 1:tm.N
		xi, yi = ind2sub(bin_size, sample(weights))
		points[n] = [em.cell_size*(xi-.5), em.cell_size*(yi-.5)]
	end

	# find a short path heuristically
	#tsp_rand!(xd, points)
	tsp_nn!(xd, points)

	# compute controls
	ud = compute_controls(xd, tm.h)

	return xd, ud
end
