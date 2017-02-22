=========================
Trajectory Manager
=========================
The `TrajectoryManager` contains information used during trajectory generation.

Dynamics
===========
Linear Dynamics are common and easy to set up:
::

    ld = LinearDynamics(A,B)

The Dubins car is a standard model for cars.
::

    tm.dynamics = DubinsDynamics(v0, r)

If you want to implement your own dynamics, you need to implement the following functions for it:
::

    function linearize(ld::DubinsDynamics, x::VF, u::VF, h::Float64)

Also,
::
    function forward_euler(dd::DubinsDynamics, x::VF, u::VF, h::Float64)


Initializer
============
::
    
    ci = ConstantInitializer(action::Vector{Float64})


Descender
============

