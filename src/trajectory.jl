######################################################################
# trajectory.jl
#
# handles the generation of ergodic trajectories
######################################################################

type TrajectoryManager

	# needed for all trajectories
	N::Int
	h::Float64
	x0::Vector{Float64}
	T::Float64
	
	# needed for ergodic trajectories
	Q::Matrix{Float64}
	R::Matrix{Float64}
	q::Float64
	max_iters::Int
	initializer::Initializer

	# dynamics stuff
	n::Int
	m::Int
	A::Matrix{Float64}
	B::Matrix{Float64}

	function TrajectoryManager(x0::Vector{Float64}, h::Float64, N::Int)
		return TrajectoryManager(x0, h, N, RandomInitializer())
	end
	function TrajectoryManager(x0::Vector{Float64}, h::Float64, N::Int, i::Initializer)
		tm = new()

		# needed for all trajectories
		tm.N = N
		tm.h = h
		tm.T = N*h
		tm.x0 = deepcopy(x0)

		# needed for ergodic trajectories
		tm.Q = eye(2)
		tm.R = 0.01 * eye(2)
		tm.q = 1.0
		tm.max_iters = 30
		tm.initializer = i

		# dynamics stuff
		tm.n = 2
		tm.m = 2
		tm.A = eye(2)
		tm.B = tm.h * eye(2)

		return tm
	end
end

"""
`dynamics!(tm::TrajectoryManager, A::Matrix{Float64}, B::Matrix{Float64})`

Sets the dynamics for `tm`.
This includes fields `tm.A`, `tm.B`, `tm.n` and `tm.m`.
Dynamics are assumed to be linear and constant.
"""
function dynamics!(tm::TrajectoryManager, A::Matrix{Float64}, B::Matrix{Float64})
	tm.n, tm.m = size(B)
	tm.A = deepcopy(A)
	tm.B = deepcopy(B)
end


function initialize(em::ErgodicManager, tm::TrajectoryManager)
	initialize(tm.initializer, em, tm)
end

# TODO: change rand() call so we can account for various
function initialize(ri::RandomInitializer, em::ErgodicManager, tm::TrajectoryManager)
	xd = [[tm.x0[1], tm.x0[2]] for i = 1:tm.N+1]
	points = Array(Vector{Float64}, tm.N)
	for i = 1:tm.N
		points[i] = [rand(), rand()]
	end

	tsp_nn!(xd, points)

	ud = compute_controls(xd, tm.h)

	return xd, ud
end



# Assumes we our domain is 2d
# corner 1 = 0,0
# corner 2 = 1,0
# corner 3 = 0,1
# corner 4 = 1,1
#function initialize(ci::CornerInitializer, em::ErgodicManager, x0::Vector{Float64}, h::Float64, N::Int)
function initialize(ci::CornerInitializer, em::ErgodicManager, tm::TrajectoryManager)

	# dimensionality of state space
	x0 = tm.x0
	n = length(x0)

	# first corner: 0,0
	dx1 = 0.0 - x0[1]
	dy1 = 0.0 - x0[2]
	dc1 = sqrt(dx1*dx1 + dy1*dy1)

	# second corner: 1,0
	dx2 = 1.0 - x0[1]
	dy2 = 0.0 - x0[2]
	dc2 = sqrt(dx2*dx2 + dy2*dy2)

	# third corner: 0,1
	dx3 = 0.0 - x0[1]
	dy3 = 1.0 - x0[2]
	dc3 = sqrt(dx3*dx3 + dy3*dy3)

	# fourth corner: 1,1
	dx4 = 1.0 - x0[1]
	dy4 = 1.0 - x0[2]
	dc4 = sqrt(dx4*dx4 + dy4*dy4)

	bi = indmax([dc1,dc2,dc3,dc4])

	bc = [0.,0.]
	if bi == 2
		bc = [1.,0.]
	elseif bi == 3
		bc = [0.,1.]
	elseif bi == 4
		bc = [1.,1.]
	end

	x_step = (bc[1] - x0[1]) / tm.N
	y_step = (bc[2] - x0[2]) / tm.N

	xd = Array(Vector{Float64}, tm.N+1)
	for i = 0:tm.N
		xd[i+1] = zeros(n)
		xd[i+1][1] = x0[1] + i*x_step
		xd[i+1][2] = x0[2] + i*y_step
	end

	ud = compute_controls(xd, tm.h)

	return xd, ud
end


#function initialize(ci::ConstantInitializer, em::ErgodicManager, x0::Vector{Float64}, h::Float64, N::Int)
function initialize(ci::ConstantInitializer, em::ErgodicManager, tm::TrajectoryManager)
	xd = Array(Vector{Float64}, tm.N+1)
	ud = Array(Vector{Float64}, tm.N+1)
	xd[1] = deepcopy(tm.x0)
	ud[1] = ci.action
	for i = 1:tm.N
		xd[i+1] = tm.A*xd[i] + tm.B*ud[i]
		ud[i+1] = deepcopy(ci.action)
	end

	return xd, ud
end



# currently, I just move in a random direction
# perhaps find direction to mean and move towards it
function initialize_trajectory(N::Int, h::Float64, x0::T2F)
	xd = [[x0[1], x0[2]] for i = 1:N+1]
	#ud = [.01*ones(2) for i = 1:N+1]
	#ud = [[.01,0.] for i = 1:N+1]
	ud = [[.001,0.001] for i = 1:N+1]
	for i = 2:N+1
		xd[i][1] = xd[i-1][1] + h*ud[i-1][1]
		xd[i][2] = xd[i-1][2] + h*ud[i-1][2]
	end
	return xd, ud
end


# computes controls from a trajectory
# TODO: really, this is a general tool useful for other code
#  it should really go somewhere else
function compute_controls(xd::VV_F, h::Float64)
	N = length(xd) - 1
	ud = Array(Vector{Float64}, N+1)
	for n = 1:N
		ud[n] = (xd[n+1] - xd[n]) / h
	end
	ud[N+1] = ud[N]   # this one doesn't matter
	
	return ud
end

# uses nearest neighbors to heuristically solve tsp
# tsp is traveling salesman problem
function tsp_nn!(xd::VV_F, points::VV_F)
	xc = xd[1]
	next_p = xd[1]
	n = 1
	while length(points) > 0
		best_d = Inf
		best_ind = 1
		for (p_ind,p) in enumerate(points)
			r = xc - p 
			d = dot(r, r)
			if d < best_d
				best_d = d
				best_ind = p_ind
			end
		end
		xd[n+1] = points[best_ind]
		deleteat!(points, best_ind)
		n += 1
	end
end

