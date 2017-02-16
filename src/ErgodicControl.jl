module ErgodicControl

# export functions I've made
export ErgodicManager, phik!, reconstruct, decompose
export ergodic_score, control_score, total_score
export control_effort
export centroid, covariance
export TrajectoryManager, dynamics!, make_trajectory, sample_trajectory
export kmeans_trajectory, random_trajectory
export assign_step, mean_step
export Initializer, initialize
export RandomInitializer, CornerInitializer, ConstantInitializer, CornerConstantInitializer, SampleInitializer, GreedyInitializer, PointInitializer
export DirectionInitializer
export clerc_trajectory
export cerc_trajectory
export controls2trajectory

# to make some things easier
typealias VV_F   Vector{Vector{Float64}}   # vector of vector of floats
typealias V_T2F  Vector{NTuple{2,Float64}} # vector of tuples of 2 floats
typealias VMF64  Vector{Matrix{Float64}}   # vector of matrix of floats
typealias VM_F  Vector{Matrix{Float64}}   # vector of matrix of floats
typealias T2F      NTuple{2, Float64}    # x, y

typealias MF Matrix{Float64}

# math-type stuff I might need
include("math.jl")
include("tsp.jl")
include("lq.jl")
include("lqr.jl")

# ergodic manager
include("ergodicity.jl")

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

end # module
