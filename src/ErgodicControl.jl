module ErgodicControl

# import functions I've made
export ErgodicManager, phik!, reconstruct, decompose, ergodic_score
export centroid, covariance
export TrajectoryManager, make_trajectory, sample_trajectory

# to make some things easier
typealias VV_F   Vector{Vector{Float64}}   # vector of vector of floats
typealias V_T2F  Vector{NTuple{2,Float64}} # vector of tuples of 2 floats
typealias VMF64  Vector{Matrix{Float64}}   # vector of matrix of floats
typealias T2F      NTuple{2, Float64}    # x, y

include("math.jl")

include("ergodicity.jl")
include("trajectory.jl")
include("max_trajectory.jl")
include("sample_trajectory.jl")

end # module
