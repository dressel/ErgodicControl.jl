######################################################################
# initializer.jl
# different ways to initalize a trajectory
######################################################################

abstract Initializer

"""
`ri = CornerInitializer()`

Takes a trajectory to the farthest corner.
"""
type RandomInitializer <: Initializer end


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
