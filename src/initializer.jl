######################################################################
# initializer.jl
# different ways to initalize a trajectory
######################################################################

abstract Initializer

"""
`ri = RandomInitializer()`

Random selection of points.
"""
type RandomInitializer <: Initializer end


"""
`si = SampleInitializer()`

Samples points from a distribution.
"""
type SampleInitializer <: Initializer end


"""
`ci = CornerInitializer()`

Takes a trajectory to the farthest corner.
"""
type CornerInitializer <: Initializer end


"""
`ci = ConstantInitializer(action::Vector{Float64})`

Just takes a constant action.
"""
type ConstantInitializer <: Initializer
	action::Vector{Float64}
end


# TODO: finish this one
"""
`cci = CornerConstantInitializer(action::Vector{Float64})`

Just takes a constant action.
"""
type CornerConstantInitializer <: Initializer
	magnitude::Float64
end

"""
`gi = GreedyInitializer()`

Greedily goes to spot with maximum phi.
Assumes phi decreases at a constant rate.
"""
type GreedyInitializer <: Initializer end


"""
`poi = PointInitializer(xd::Vector{Float64}, mag::Float64)`

Moves in the direction of point `xd` with magnitude `mag`.
"""
type PointInitializer <: Initializer
	xd::Vector{Float64}
	mag::Float64
end
