=========================
Trajectory Manager
=========================
The :code:`TrajectoryManager` contains information used during trajectory generation.

Fields and Construction
=========================
The :code:`TrajectoryManager` type has the following fields:
::

	# basic info
	N::Int                      # trajectory length
	h::Float64                  # time step
	x0::Vector{Float64}         # initial state
	
	# cost functions
	q::Float64                  # ergodic cost multiplier, default = 1
	R::Matrix{Float64}          # control cost multiplier, def = .01*eye(2)
	Qn::Matrix{Float64}         # LQ ergodic cost, default = eye(2)
	Rn::Matrix{Float64}         # LQ control cost, default = eye(2)
	barrier_cost::Float64       # penalizes leaving domain, def = 0

	# needed for trajectory generation
	initializer::Initializer
	descender::Descender        # default is ArmijoLineSearch(10,.1)
	dynamics::Dynamics          # default is single integrator

A :code:`TrajectoryManager` is constructed with basic information about the trajectory and an optional initializer.
::

	TrajectoryManager(x0::Vector{Float64}, h::Float64, N::Int, i::Initializer=RandomInitializer())

By Default,  :code:`barrier_cost=0`, meaning no barrier cost is applied. When :code:`barrier_cost` is positive, a quadratic barrier function is added to the objective. This cost penalizes the trajectory for leaving the domain specified in the :code:`ErgodicManager` during trajectory generation.

The :code:`initializer`, :code:`descender`, and :code:`dynamics` fields are described below.


Initializer
============
The :code:`initializer` field must be a subtype of the abstract `Initializer` type.

A good option is the :code:`ConstantInitializer`, which just sets every action to the provided value and uses a forward Euler method to determine the trajectory. Below is an example when the control inputs are in :math:`\mathbb{R}^2`.
::
    
    tm.initializer = ConstantInitializer([0.0,0.0])


Descender
============
The :code:`descender` field must be a subtype of the :code:`Descender` abstract type.

The default is an :code:`ArmijoLineSearch`, as Armijo line search has been used extensively in the literature.

Simpler alternatives are shown below. They work ok, but seem to take longer than Armijo line search (which is to be expected). I'll use them when I'm troubleshooting or if Armijo line search doesn't work well in a particular problem.
::

    tm.descender = ConstantStep(0.5)        # each step is 0.5
    tm.descender = InverseStep(0.5)         # each step is 0.5/iter
    tm.descender = InverseRootStep(0.15)    # each step is 0.15/sqrt(iter)


Dynamics
===========
Linear dynamics are common and easy to set up:
::

    tm.dynamics = LinearDynamics(A,B)

The Dubins car is a standard model for cars.
::

    tm.dynamics = DubinsDynamics(v0, r)

If you want to implement your own dynamics, you need to subtype the abstract :code:`Dynamics` type and implement the :code:`linearize` and :code:`forward_euler` functions. 
::

    type MyDynamics <: Dynamics
        # fields, constructors, etc

        # following fields must be included
        n::Int      # number of state variables
        m::Int      # number of control variables
    end

    function linearize(md::MyDynamics, x::VF, u::VF, h::Float64)
        # return A,B matrices
    end

    function forward_euler(md::MyDynamics, x::VF, u::VF, h::Float64)
        # return new state
    end
