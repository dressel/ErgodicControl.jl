######################################################################
# trajectory_manager.jl
#
# handles the trajectory manager
######################################################################

using StatsBase: WeightVec, sample

abstract Dynamics
abstract Initializer
abstract Descender

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

	# The descender determines how much to descend at each step
	descender::Descender

	# dynamics stuff
	dynamics::Dynamics

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
		# trying a higher cost now 2/07/2017
		#tm.R = 0.01 * eye(2)
		tm.R = 0.01 * eye(2)
		tm.q = 1.0
		tm.max_iters = 30
		tm.initializer = i
		tm.descender = InverseRootStep(1.0)
		# I've found InverseRootStep(0.15) to work well

		# dynamics stuff
		tm.dynamics = LinearDynamics(eye(2), tm.h*eye(2))

		return tm
	end
end


"""
`dynamics!(tm::TrajectoryManager, A::Matrix{Float64}, B::Matrix{Float64})`

`dynamics!(tm::TrajectoryManager, d_type::String)`

Sets the dynamics for `tm`.
This includes fields `tm.A`, `tm.B`, `tm.n` and `tm.m`.
Dynamics are assumed to be linear and constant.

User can optionally pass in string `d_type`.
If this string is "double integrator", double integrator dynamics will be made.

If the new `tm.n` (number of state dimensions) is greater than the length of `tm.x0`, `tm.x0` will be padded with zeros until it is of length `tm.n`.

"""
function dynamics!(tm::TrajectoryManager, A::Matrix{Float64}, B::Matrix{Float64})
	ld = LinearDynamics(A, B)
	tm.dynamics = ld
	while length(tm.x0) < ld.n
		push!(tm.x0, 0.)
	end
end

function dynamics!(tm::TrajectoryManager, d_type::String)
	if d_type == "double integrator"
		A = [1 0 tm.h 0; 0 1 0 tm.h; 0 0 1 0; 0 0 0 1]
		B = [0 0; 0 0; tm.h 0; 0 tm.h]
		dynamics!(tm, A, B)
	elseif d_type == "dubins"
		tm.dynamics = DubinsDynamics(1.0, 1.0)
	else
		error("Invalid d_type provided!")
	end
end


function initialize(em::ErgodicManager, tm::TrajectoryManager)
	initialize(tm.initializer, em, tm)
end

# TODO: change rand() call so we can account for various
export compute_controls


# computes controls from a trajectory
# TODO: really, this is a general tool useful for other code
#  it should go somewhere else
function compute_controls(xd::VV_F, h::Float64)
	N = length(xd) - 1
	ud = Array(Vector{Float64}, N)
	for n = 1:N
		ud[n] = (xd[n+1] - xd[n]) / h
	end
	
	return ud
end

"""
controls2trajectory(tm::TrajectoryManager, ud::VV_F)

computes trajectory from controls (and initial position)
"""
function controls2trajectory(tm::TrajectoryManager, ud::VV_F)
	N = length(ud)
	xd = Array(Vector{Float64}, N+1)

	xd[1] = deepcopy(tm.x0)
	for i = 1:tm.N
		xd[i+1] = tm.dynamics.A*xd[i] + tm.dynamics.B*ud[i]
	end

	return xd
end


# creates a sample_trajectory
function sample_trajectory(em::ErgodicManager, tm::TrajectoryManager)
	return initialize(SampleInitializer(), em, tm)
end



function linearize(tm::TrajectoryManager, x::Vector{Float64}, u::Vector{Float64})
	return linearize(tm.dynamics, x, u, tm.h)
end
function forward_euler(tm::TrajectoryManager, x::Vector{Float64}, u::Vector{Float64})
	forward_euler(tm.dynamics, x, u, tm.h)
end
