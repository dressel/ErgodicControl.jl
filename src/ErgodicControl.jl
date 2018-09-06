__precompile__()

module ErgodicControl

import Base.normalize!
using SpecialFunctions

# export functions I've made
export ErgodicManager, decompose!, reconstruct, decompose
export ergodic_score, control_score, total_score
export control_effort
export centroid, covariance
export TrajectoryManager, dynamics!, sample_trajectory
export kmeans_trajectory, random_trajectory
export assign_step, mean_step

export
    Initializer,
    initialize,
    CornerConstantInitializer,
    GreedyInitializer,
    PointInitializer,
    DirectionInitializer

export
    x_min, y_min, z_min,
    x_max, y_max, z_max,
    x_size, y_size, z_size,
    x_cells, y_cells, z_cells

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

end # module
