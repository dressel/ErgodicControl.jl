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
	
	# Cost functions
	q::Float64
	Qn::Matrix{Float64}
	R::Matrix{Float64}
	Rn::Matrix{Float64}
	barrier_cost::Float64

	initializer::Initializer
	descender::Descender
	dynamics::Dynamics

	#function TrajectoryManager(x0::Vector{Float64}, h::Float64, N::Int)
	#	return TrajectoryManager(x0, h, N, RandomInitializer())
	#end
	function TrajectoryManager(x0::Vector{Float64}, h::Real, N::Int, i::Initializer=RandomInitializer())
		tm = new()

		# needed for all trajectories
		tm.N = N
		tm.h = h
		tm.T = N*h
		tm.x0 = deepcopy(x0)

		# needed for ergodic trajectories
		tm.Qn = eye(2)
		tm.q = 1.0
		tm.R = 0.01 * eye(2)
		tm.Rn = eye(2)
		tm.barrier_cost = 0.
		tm.initializer = i
		#tm.descender = InverseRootStep(1.0)
		tm.descender = ArmijoLineSearch()

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



# computes controls from a trajectory
# TODO: really, this is a general tool useful for other code
#  it should go somewhere else
export compute_controls
function compute_controls(xd::VV_F, h::Float64)
	N = length(xd) - 1
	ud = Array(Vector{Float64}, N)
	for n = 1:N
		ud[n] = (xd[n+1] - xd[n]) / h
	end
	
	return ud
end


# creates a sample_trajectory
function sample_trajectory(em::ErgodicManager, tm::TrajectoryManager)
	return initialize(SampleInitializer(), em, tm)
end
