=========================
Examples 
=========================

Single Integrator
==================
Here is an example
::

    using ErgodicControl

    em = ErgodicManagerR2("single gaussian", K=5, bins=100)

    x0 = [0.4,0.1]
    N = 40
    h = 0.1

    tm = TrajectoryManager(x0, h, N, ConstantInitializer([0.0,0.0]))

    xd, ud = pto_trajectory(em, tm)
    plot(em, xd)

.. image:: http://stanford.edu/~dressel/gifs/ergodic/single1.png



Double Integrator
===================
Here is another example
::

    using ErgodicControl

    em = ErgodicManagerR2("double gaussian", K=5, bins=100)

    x0 = [0.49,0.01,0.0,0.0]
    N = 40
    h = 0.5

    tm = TrajectoryManager(x0, h, N, ConstantInitializer([0.0,0.0]))

    # dynamics stuff
    dynamics!(tm, "double integrator")
    tm.Qn = eye(4)

    tm.descender = ArmijoLineSearch(1,.01)

    xd, ud = pto_trajectory(em, tm)

When doing this, the solver gets stuck. We can try another descent engine.
::

    tm.descender = ArmijoLineSearch(1,.01)

.. image:: http://stanford.edu/~dressel/gifs/ergodic/double1.png


Dubins Car
===================
The Dubins car is simple and often used to model cars.
::

    using ErgodicControl

    em = ErgodicManagerR2("double gaussian", K=5, bins=100)

    x0 = [0.5,0.01,pi/4]
    N = 40
    h = 0.1

    tm = TrajectoryManager(x0, h, N, ConstantInitializer([0.0000]))

    # things needed for dynamics
    tm.dynamics = DubinsDynamics(0.3,0.1)
    tm.Qn = eye(3)
    tm.R = 0.01 * eye(1)
    tm.Rn = 1 * eye(1)

    xd, ud = pto_trajectory(em, tm)
    plot(em, xd)

Look at how bad this is! The vehicle leaves our domain!

.. image:: http://stanford.edu/~dressel/gifs/ergodic/dubins1.png

As mentioned before, the trajectory leaves the domain because the Fourier basis function is periodic. This makes sense in the context of the Dubins dynamics. Control effort is only expended when changing the vehicle's heading.

We can overcome this by penalizing states outside the domain, using the barrier cost we mentioned before. This is handled with the trajectory manager's barrier_cost field, which is set to 0 by default. Let's try changing the cost to 1. We'd add the following line before the call to :code:`pto_trajectory`:
::

    tm.barrier_cost = 1

With this modification, trajectory generation reaches the directional derivative criterion after 155 iterations. The trajectory stays within the domain.

.. image:: http://stanford.edu/~dressel/gifs/ergodic/dubins2.png


Distributions in Three Dimensions
==================================


Time-evolving Spatial Distribution
========================================


Multi-agent Trajectories
===============================
This works
::

    using ErgodicControl

    # Set up different domains with different discretizations
    d = Domain([1,1], 100)
    num_agents = 2

    # Set up distribution and ergodic manager
    K = 5
    means = [[.3,.7], [.7,.3]]
    Sigmas = [.025*eye(2), .025*eye(2)]
    phi = gaussian(d, means, Sigmas)
    em = ErgodicManagerR2(d, phi, K)

    # Set up first trajectory manager
    x0 = [0.49,0.01]
    N = 50
    h = 0.6
    ci = ConstantInitializer([0.0, 0.0])
    tm1 = TrajectoryManager(x0, h, N, ci)
    dynamics!(tm1, SingleIntegrator(2,h))

    # second tm is like the first, but different starting point
    tm2 = deepcopy(tm1)
    tm2.x0 = [.79,.99]

    # array of trajectory managers
    vtm = [tm1, tm2]

    # Generate the trajectories
    ddc = 1e-4
    xd, ud = pto_trajectory(em, vtm, dd_crit=ddc)
    plot(em, xd, vtm)


.. image:: http://stanford.edu/~dressel/gifs/ergodic/multi.png


Distribution over SE(2)
===============================
