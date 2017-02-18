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


function initialize(ci::ConstantInitializer, em::ErgodicManager, tm::TrajectoryManager)
	xd = Array(Vector{Float64}, tm.N+1)
	ud = Array(Vector{Float64}, tm.N)
	xd[1] = deepcopy(tm.x0)
	ud[1] = deepcopy(ci.action)
	for i = 1:(tm.N-1)
		# 2/07/2017: need to compute true A
		#A,B = linearize(tm, xd[i], ud[i])
		#display(A)
		#xd[i+1] = A*xd[i] + B*ud[i]
		xd[i+1] = forward_euler(tm, xd[i], ud[i])
		ud[i+1] = deepcopy(ci.action)
	end
	#A,B = linearize(tm, xd[tm.N], ud[tm.N])
	#xd[tm.N+1] = A*xd[tm.N] + B*ud[tm.N]
	xd[tm.N+1] = forward_euler(tm, xd[tm.N], ud[tm.N])

	return xd, ud
end
