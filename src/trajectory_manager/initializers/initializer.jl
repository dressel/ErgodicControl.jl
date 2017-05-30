######################################################################
# initializer.jl
# different ways to initalize a trajectory
######################################################################
export ConstantInitializer
export RandomInitializer
export SampleInitializer
export CornerInitializer

include("constant_initializer.jl")
include("random_initializer.jl")
include("sample_initializer.jl")
include("corner_initializer.jl")


"""
`xd, ud = initialize(em::ErgodicManager, tm::TrajectoryManager)`

Runs `tm`'s initializer to return a trajectory.
"""
function initialize(em::ErgodicManager, tm::TrajectoryManager)
	initialize(tm.initializer, em, tm)
end
function initialize(em::ErgodicManager, vtm::Vector{TrajectoryManager})
	num_agents = length(vtm)
	xd0s = VVVF()
	ud0s = VVVF()
	for j = 1:num_agents
		xd0, ud0 = initialize(em, vtm[j])
		push!(xd0s, xd0)
		push!(ud0s, ud0)
	end
	return vvvf2vvf(xd0s, ud0s)
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
`poi = PointInitializer(xd::Vector{Float64})`

Moves to point `xd`.
"""
type PointInitializer <: Initializer
	xd::Vector{Float64}
end

"""
`di = DirectionInitializer(xd::Vector{Float64}, mag::Float64)`

Moves in the direction of point `xd` with magnitude `mag`.
"""
type DirectionInitializer <: Initializer
	xd::Vector{Float64}
	mag::Float64
end





# TODO: actually finish this
function initialize(cci::CornerConstantInitializer, em::ErgodicManager, tm::TrajectoryManager)
	xd = Array(Vector{Float64}, tm.N+1)
	ud = Array(Vector{Float64}, tm.N)
	xd[1] = deepcopy(tm.x0)
	ud[1] = deepcopy(ci.action)
	for i = 1:(tm.N-1)
		xd[i+1] = tm.A*xd[i] + tm.B*ud[i]
		ud[i+1] = deepcopy(ci.action)
	end
	xd[tm.N+1] = tm.A*xd[tm.N] + tm.B*ud[tm.N]

	return xd, ud
end


export greedy_trajectory
function greedy_trajectory(em::ErgodicManager, tm::TrajectoryManager)
	return initialize(GreedyInitializer(), em, tm)
end
function initialize(gi::GreedyInitializer, em::ErgodicManager, tm::TrajectoryManager)
	d_rate = sum(em.phi)/tm.N
	num_cells = em.bins*em.bins
	total_info = 0.0
	xd = Array(Vector{Float64}, tm.N+1)
	xd[1] = deepcopy(tm.x0)
	temp_phi = deepcopy(em.phi)
	size_tuple = (em.bins, em.bins)
	for n = 1:tm.N
		bi = indmax(temp_phi)
		xi, yi = ind2sub(size_tuple, bi)
		xd[n+1] = [(xi-0.5)*em.cell_size, (yi-0.5)*em.cell_size]
		temp_phi[bi] -= min(temp_phi[bi], d_rate)
	end
	ud = compute_controls(xd, tm.h)
	return xd,ud
end


# moves to a point with a constant control input
function initialize(initializer::PointInitializer, em::ErgodicManager, tm::TrajectoryManager)

	# compute the direction and action we most go towards
	dx = initializer.xd[1] - tm.x0[1]
	dy = initializer.xd[2] - tm.x0[2]
	#x_step = (initializer.xd[1] - tm.x0[1]) / (tm.N * tm.h)
	#y_step = (initializer.xd[2] - tm.x0[2]) / (tm.N * tm.h)
	#u = [x_step, y_step]

	# if we are double integrator, only apply the input once
	ud = Array(Vector{Float64}, tm.N)
	if tm.dynamics.n > 2
		den = tm.h * tm.h * (tm.N-1)
		ud[1] = [dx/den, dy/den]
		for i = 1:(tm.N-1)
			ud[i+1] = zeros(2)
		end
	else  # single integrator
		den = tm.h * tm.N
		u = [dx / den, dy / den]
		for i = 1:tm.N
			ud[i] = deepcopy(u)
		end
	end

	xd = integrate(tm, ud)

	return xd, ud
end

# moves to a point with a constant control input
function initialize(initializer::DirectionInitializer, em::ErgodicManager, tm::TrajectoryManager)
	xd = Array(Vector{Float64}, tm.N+1)
	ud = Array(Vector{Float64}, tm.N)

	# compute the direction and action we most go towards
	dx = initializer.xd[1] - tm.x0[1]
	dy = initializer.xd[2] - tm.x0[2]
	u = initializer.mag * [dx, dy] / sqrt(dx*dx + dy*dy)

	xd[1] = deepcopy(tm.x0)
	ud[1] = deepcopy(u)
	for i = 1:(tm.N-1)
		xd[i+1] = tm.A*xd[i] + tm.B*ud[i]
		ud[i+1] = deepcopy(u)
	end
	xd[tm.N+1] = tm.A*xd[tm.N] + tm.B*ud[tm.N]

	return xd, ud
end
