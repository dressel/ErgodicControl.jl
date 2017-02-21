=========================
Examples 
=========================

Single Integrator
==================
Here is an example
::

    using ErgodicControl

    em = ErgodicManager("single gaussian", K=5, bins=100)

    x0 = [0.4,0.1]
    N = 40
    h = 0.1

    tm = TrajectoryManager(x0, h, N, ConstantInitializer([0.0,0.0]))

    @time xd, ud = new_trajectory(em, tm, logging=true)
    #gif(em, tm, fps=15)
    plot(em, xd)

.. image:: http://stanford.edu/~dressel/gifs/gif3.gif


Double Integrator
===================
Here is another example
::
    
    ok there


Dubins Car
===================
The Dubins car is blah blah
::

    using ErgodicControl

    em = ErgodicManager("double gaussian", K=5, bins=100)

    x0 = [0.5,0.01,pi/4]
    N = 40
    h = 0.1

    tm = TrajectoryManager(x0, h, N, ConstantInitializer([0.0000]))

    # things needed for dynamics
    tm.dynamics = DubinsDynamics(0.3,0.1)
    tm.Qn = eye(3)
    tm.R = 0.01 * eye(1)
    tm.Rn = 1 * eye(1)

    @time xd, ud = new_trajectory(em, tm, max_iters=1000)
    plot(em, xd)

One way we can handle this is
::

    tm.barrier_cost = 1
