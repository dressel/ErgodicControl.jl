# using ErgodicControl
# using Test
#
# @testset "ErgodicControl.jl" begin
#     # Write your tests here.
# end

using ErgodicControl
using Base.Test

# write your own tests here
@test 1 == 1

em = ErgodicManagerR2("single gaussian", K=5, bins=100)
x0 = [0.4,0.1]
N = 40
h = 0.1

tm = TrajectoryManager(x0, h, N, ConstantInitializer([0.0,0.0]))
xd, ud = pto_trajectory(em, tm)

include("dubins.jl")

include("r2t.jl")

include("multi.jl")

include("r3.jl")

include("smc.jl")

include("se2.jl")

include("linear.jl")

include("max.jl")
