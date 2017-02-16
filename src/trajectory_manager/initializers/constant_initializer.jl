######################################################################
# constant_initialier.jl
######################################################################
"""
`ci = ConstantInitializer(action::Vector{Float64})`

Just takes a constant action.
"""
type ConstantInitializer <: Initializer
	action::Vector{Float64}
end
