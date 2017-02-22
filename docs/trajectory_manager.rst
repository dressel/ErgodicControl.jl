=========================
Trajectory Manager
=========================
The `TrajectoryManager` contains information used during trajectory generation.

Fields and Construction
=========================
The `TrajectoryManager` type has the following fields:
::

	# needed for all trajectories
	N::Int                      # trajectory length
	h::Float64                  # time step
	x0::Vector{Float64}         # initial state
	T::Float64                  # time horizon (N * h)
	
	# Cost functions
	q::Float64                  # ergodic score multiplier
	Qn::Matrix{Float64}         # LQ ergodic cost 
	R::Matrix{Float64}          # control score multiplier
	Rn::Matrix{Float64}         # LQ control cost
	barrier_cost::Float64

	initializer::Initializer
	descender::Descender
	dynamics::Dynamics

A `TrajectoryManager` is constructed with basic information about the trajectory and an optional initializer.
::

	TrajectoryManager(x0::Vector{Float64}, h::Float64, N::Int, i::Initializer=RandomInitializer())

The default

The `barrier_cost` field is set to 0 by default, meaning no barrier cost is applied. When `barrier_cost` is positive, a quadratic barrier function is added to the objective.

The `initializer`, `descedner`, and `dynamics` fields are described below.


Initializer
============
::
    
    ci = ConstantInitializer(action::Vector{Float64})


Descender
============
The `descender`


Dynamics
===========
Linear Dynamics are common and easy to set up:
::

    tm.dynamics = LinearDynamics(A,B)

The Dubins car is a standard model for cars.
::

    tm.dynamics = DubinsDynamics(v0, r)

If you want to implement your own dynamics, you need to subtype the abstract `Dynamics` type and implement the `linearize` and `forward_euler` functions.
::

    type MyDynamics <: Dynamics
        # fields, constructors, etc
    end

    function linearize(md::MyDynamics, x::VF, u::VF, h::Float64)
        # return A,B matrices
    end

    function forward_euler(md::MyDynamics, x::VF, u::VF, h::Float64)
        # return new state
    end


