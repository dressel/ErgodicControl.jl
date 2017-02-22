==========================
Generating Trajectories
==========================

Once you have created an ErgodicManager and a TrajectoryManager, you can generate a trajectory.
::

    xd, ud = new_trajectory(em, tm)

There are a number of optional arguments. Here they are with their default arguments:
::

    verbose::Bool = true
    logging::Bool = false
    max_iters::Int = 100
    es_crit::Float64 = 0.003
    dd_crit::Float64 = 1e-6
    right::Bool = false

When the `verbose` tag is set to `true`, progress is presented at each descent iteration. You can turn this to false if you are running trajectory generation as an inner component of a larger algorithm.

The `logging` tag saves a file `temp.csv` that contains a copy of each trajectory (just the states, not the actions) at each descent iteration. This is mostly there for gif generation (see Visuals).

The `max_iters` tag is the maximum number of descent iterations allowed.

The `es_crit` and `dd_crit` tags are termination conditions based on the ergodic score and directional derivative, respectively. A lot of research suggests using a directional derivative criterion.

The `right` tag determines if right or left-Riemann sums are used when generating trajectories. Right now, only the `false` option works (that is, using left-Riemann sums).
