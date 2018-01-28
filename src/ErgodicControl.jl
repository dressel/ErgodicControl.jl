module ErgodicControl

import Base.normalize!

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
# deprecation
#typealias T2F      NTuple{2, Float64}    # x, y
#typealias MF  Matrix{Float64}
#typealias VMF Vector{MF}
#typealias VF  Vector{Float64}
#typealias VVF Vector{VF}
#typealias VVVF Vector{VVF}

# to make some things easier
const T2F = NTuple{2, Float64}    # x, y
const MF = Matrix{Float64}
const VMF = Vector{MF}
const VF = Vector{Float64}
const VVF = Vector{VF}
const VVVF = Vector{VVF}

# math-type stuff I might need
include("math.jl")
include("tsp.jl")
include("lq.jl")
include("lqr.jl")

# ergodic manager
include("ergodic_manager/domain.jl")
include("ergodic_manager/gaussian.jl")
include("ergodic_manager/circle.jl")
include("ergodic_manager/ergodic_manager.jl")
include("ergodic_manager/r2.jl")
include("ergodic_manager/r3.jl")
include("ergodic_manager/r2t.jl")
include("ergodic_manager/se2.jl")
include("ergodic_manager/examples.jl")

# trajectory manager stuff
include("trajectory_manager/trajectory_manager.jl")
include("trajectory_manager/dynamics.jl")
include("trajectory_manager/initializers/initializer.jl")
include("trajectory_manager/descender.jl")

# trajectory generation
include("trajectory_generation/trajectory_generation.jl")

# include visuals
#include("visuals/plots.jl")
#include("visuals/gif.jl")


end # module
