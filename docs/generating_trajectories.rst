==========================
Generating Trajectories
==========================

Once you have created an ErgodicManager and a TrajectoryManager, you can generate a trajectory.
::

    cerc_trajectory(em, tm)

There are a number of optional arguments:
::

    verbose::Bool = true
    logging::Bool = false
    max_iters::Int = 30
    right::Bool =true
