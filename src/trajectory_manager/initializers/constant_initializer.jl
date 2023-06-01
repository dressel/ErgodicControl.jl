######################################################################
# constant_initialier.jl
######################################################################
"""
`ci = ConstantInitializer(action::Vector{Float64})`

Just takes a constant action.
"""
mutable struct ConstantInitializer <: Initializer
	action::Vector{Float64}
end


function initialize(ci::ConstantInitializer, em::ErgodicManager, tm::TrajectoryManager)
	xd =  Array{Vector{Float64}}(undef, tm.N+1)#1.0 .*ones(tm.N+1,2,2)#Array{Vector{Float64}}(tm.N+1)
	ud =  Array{Vector{Float64}}(undef, tm.N)#1.0 .*ones(tm.N,2,2)#Array{Vector{Float64}}(tm.N)
	xd[1] = deepcopy(tm.x0)
	ud[1] = deepcopy(ci.action)
	for i = 1:(tm.N-1)
		xd[i+1] = integrate(tm, xd[i], ud[i])
		ud[i+1] = deepcopy(ci.action)
	end
	xd[tm.N+1] = integrate(tm, xd[tm.N], ud[tm.N])

	return xd, ud
end
