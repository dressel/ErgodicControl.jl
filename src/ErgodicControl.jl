module ErgodicControl

# export functions I've made
export ErgodicManager, decompose!, reconstruct, decompose
export ergodic_score, control_score, total_score
export control_effort
export centroid, covariance
export TrajectoryManager, dynamics!, make_trajectory, sample_trajectory
export kmeans_trajectory, random_trajectory
export assign_step, mean_step
export Initializer, initialize
export CornerConstantInitializer, GreedyInitializer, PointInitializer
export DirectionInitializer
export clerc_trajectory
export cerc_trajectory

# to make some things easier
typealias VV_F   Vector{Vector{Float64}}   # vector of vector of floats
typealias V_T2F  Vector{NTuple{2,Float64}} # vector of tuples of 2 floats
typealias VMF64  Vector{Matrix{Float64}}   # vector of matrix of floats
typealias VM_F  Vector{Matrix{Float64}}   # vector of matrix of floats
typealias T2F      NTuple{2, Float64}    # x, y

typealias MF Matrix{Float64}
typealias VMF Vector{MF}
typealias VF   Vector{Float64}   # vector of vector of floats
typealias VVF   Vector{VF}   # vector of vector of floats

# math-type stuff I might need
include("math.jl")
include("tsp.jl")
include("lq.jl")
include("lqr.jl")

# ergodic manager
include("domain.jl")
include("ergodic_manager.jl")
include("ergodic_manager_r2.jl")
include("ergodic_manager_se2.jl")

# trajectory manager stuff
include("trajectory_manager/trajectory_manager.jl")
include("trajectory_manager/dynamics.jl")
include("trajectory_manager/initializers/initializer.jl")
include("trajectory_manager/descender.jl")

# trajectory generation
include("trajectory_generation/trajectory_generation.jl")

# include visuals
include("visuals/plots.jl")
include("visuals/gif.jl")

# temporary...
function fletcher(N::Int)
	xd = VVF(N)
	s = sqrt(0.03)
	for i = 1:N
		xd[i] = [s*randn() + 0.5, s*randn() + .5, s*randn()]
	end
	return xd
end
function nagumo(N::Int)
	xd = VVF(N)
	s = sqrt(0.03)
	for i = 1:N
		xd[i] = [rand(), randn(), rand()]
	end
	return xd
end
export nagumo
function fletcher2(N::Int, phi)
	xd = VVF(N)
	bins = size(phi,1)
	bin_size = (bins, bins)
	cell_size = 1.0 / bins

	# first sample N points from e.phi
	temp_phi = reshape(sum(phi,3), bins, bins)
	weights = WeightVec(vec(temp_phi))
	for n = 1:N
		xi, yi = ind2sub(bin_size, sample(weights))
		xd[n] = [cell_size*(xi-.5), cell_size*(yi-.5)]
	end

	return xd
end
export fletcher, fletcher2

end # module
