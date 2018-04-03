=========================
Examples 
=========================

Single Integrator
==================

::

    using ErgodicControl

    em = ErgodicManagerR2("single gaussian", K=5, bins=100)

    x0 = [0.4,0.1]
    N = 40
    h = 0.1

    tm = TrajectoryManager(x0, h, N, ConstantInitializer([0.0,0.0]))

    xd, ud = pto_trajectory(em, tm)

    # plotting
    using ErgodicControlPlots
    plot(em, xd)

.. image:: http://stanford.edu/~dressel/gifs/ergodic/example1.png



Double Integrator
===================
This example needs to be redone.
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

    # ergodic manager
    em = ErgodicManagerR2("double gaussian", K=5, bins=100)

    # trajectory manager
    x0 = [0.5,0.01,pi/4]
    N = 40
    h = 0.1
    tm = TrajectoryManager(x0, h, N, ConstantInitializer([0.0000]))

    # things needed for dynamics
    tm.dynamics = DubinsDynamics(0.3,0.1)
    tm.Qn = eye(3)
    tm.R = 0.01 * eye(1)
    tm.Rn = 1 * eye(1)

    # generate trajectory
    xd, ud = pto_trajectory(em, tm)

    # plotting
    using ErgodicControlPlots
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
An example if we have a distribution over three dimensions. The agent is a single integrator:
::

    using ErgodicControl

    # domain, distribution, and ergodic manager
    d = Domain([1,1,1], 100)
    means = [[.2,.2,.2], [.8,.8,.2], [.5,.5,.8]]
    covs = [0.01*eye(3), 0.01*eye(3), .01*eye(3)]
    phi = gaussian(d, means, covs)
    K = 5
    em = ErgodicManagerR3(d, phi, K)

    # trajectory params
    x0 = [0.49, 0.01, 0.01]
    dt = 0.5
    N = 80
    tm = TrajectoryManager(x0, dt, N, ConstantInitializer([0.0,0.0,0.0]))
    dynamics!(tm, SingleIntegrator(3,dt))
    tm.descender = ArmijoLineSearch(1,1e-4)

    # trajectory generation and plotting
    xd,ud = pto_trajectory(em, tm, dd_crit=1e-4, max_iters=1000)
    plot(em, xd, show_score=false)

Plotting is a bit trickier, and is not finished for three dimensions. The tough part is plotting the distribution. Ideally, you'd just plot some isosurfaces for the distribution, but Matplotlib wasn't made to do such things. I could try Mayavi, but that sounds like a pain. In the following image, I used a scatter plot with points sample from the distribution as a rough representation of the distribution.

.. image:: http://stanford.edu/~dressel/gifs/ergodic/three.png


Time-evolving Spatial Distribution
========================================
::

    using ErgodicControl

    # Generate the distribution
    N = 80
    dt = 0.5
    T = N*dt
    d = Domain([1,1], [100,100])
    cov = 0.010 * eye(2)
    phi = zeros(100,100,N+1)
    for i = 1:N+1
        mui = (.7*(i-1)/N + .15) * ones(2)
        phi[:,:,i] = gaussian(d, mui, cov)
    end
    ErgodicControl.normalize!(phi, d.cell_size / (N+1))

    # Now let's create the ergodic manager in R3
    K = 5
    em = ErgodicManagerR2T(d, phi, K)

    # trajectory params
    x0 = [0.49, 0.01]
    tm = TrajectoryManager(x0, dt, N, ConstantInitializer([0.,0.]))
    tm.R = .1*eye(2)

    # I call this second Armijo
    tm.descender = ArmijoLineSearch(1,1e-4)

    # trajectory generation and plotting
    mi = 1000
    ddc = 1e-5
    v = true
    xd,ud = pto_trajectory(em, tm, dd_crit=ddc, max_iters=mi, verbose=v)

    # generating the gif
    using ErgodicControlPlots
    gif(em, xd)

.. image:: http://stanford.edu/~dressel/gifs/ergodic/time.gif


Multi-agent Trajectories
===============================
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

    # plotting
    using ErgodicControlPlots
    plot(em, xd, vtm)


.. image:: http://stanford.edu/~dressel/gifs/ergodic/multi.png

Multi-agent Trajectory for Time-evolving Distribution
========================================================
We can generate a multi-agent trajectory for a time-evolving distribution.
::

    using ErgodicControl

    # Generate the distribution
    N = 80
    dt = 0.5
    T = N*dt
    d = Domain([1,1], [100,100])
    cov = 0.020 * eye(2)
    phi = zeros(100,100,N+1)
    for i = 1:N+1
        mui = (.7*(i-1)/N + .15) * ones(2)
        phi[:,:,i] = gaussian(d, mui, cov)
    end
    ErgodicControl.normalize!(phi, d.cell_size / (N+1))

    # Now let's create the ergodic manager in R2T
    K = 5
    em = ErgodicManagerR2T(d, phi, K)

    # trajectory params
    x0 = [0.49, 0.01, 0., 0.]
    tm = TrajectoryManager(x0, dt, N, ConstantInitializer([0.,0.]))
    tm.R = .01*eye(2)
    tm.descender = ArmijoLineSearch(1,1e-4)
    dynamics!(tm, DoubleIntegrator(2,dt))

    # create a vector of trajectory managers
    tm2 = deepcopy(tm)
    tm2.x0 = [.3,.9, 0., 0.]
    vtm = [tm, tm2]

    # trajectory generation and plotting
    mi = 1000
    ddc = 1e-5
    v = true
    xd,ud = pto_trajectory(em, vtm, dd_crit=ddc, max_iters=mi, verbose=v)
    gif(em, xd, vtm)

The resulting gif is shown below:

.. image:: http://stanford.edu/~dressel/gifs/ergodic/multitime.gif

The following example is also cool. The multi-agent system consists of a Dubins vehcile and a double integrator.
::

    using ErgodicControl

    # Generate the distribution
    N = 80
    dt = 0.5
    T = N*dt
    d = Domain([1,1], [100,100])
    cov = 0.020 * eye(2)
    phi = zeros(100,100,N+1)
    for i = 1:N+1
        mui = (.7*(i-1)/N + .15) * ones(2)
        phi[:,:,i] = gaussian(d, mui, cov)
    end
    ErgodicControl.normalize!(phi, d.cell_size / (N+1))

    # Now let's create the ergodic manager
    K = 5
    em = ErgodicManagerR2T(d, phi, K)

    # trajectory params
    x0 = [0.5, 0.9, 0., 0.]
    tm1 = TrajectoryManager(x0, dt, N, ConstantInitializer([0.,0.]))
    tm1.R = .01*eye(2)
    dynamics!(tm1, DoubleIntegrator(2,dt))
    tm1.barrier_cost = 1.

    tm2 = deepcopy(tm1)
    dynamics!(tm2, DubinsDynamics(.05, .1))
    tm2.initializer = ConstantInitializer([0.05])
    tm2.x0 = [.1,.1, .0]

    vtm = [tm1, tm2]

    # trajectory generation and plotting
    xd,ud = pto_trajectory(em, vtm, dd_crit=1e-5, max_iters=1000)
    gif(em, xd, vtm)

The resulting gif is shown below

.. image:: http://stanford.edu/~dressel/gifs/ergodic/dubins_doubleintegrator.gif


Distribution over SE(2)
===============================


