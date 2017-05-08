==========================
Generating Trajectories
==========================

Once you have created an ErgodicManager and a TrajectoryManager, you can generate a trajectory.
::

    xd, ud = pto_trajectory(em, tm)

There are a number of optional keyword arguments. Here they are with their default arguments:
::

    verbose::Bool = true
    logging::Bool = false
    max_iters::Int = 100
    es_crit::Float64 = 0.003
    dd_crit::Float64 = 1e-6

When the :code:`verbose` tag is set to :code:`true`, progress is presented at each descent iteration. You can turn this to false if you are running trajectory generation as an inner component of a larger algorithm.

The :code:`logging` tag saves a file :code:`temp.csv` that contains a copy of each trajectory (just the states, not the actions) at each descent iteration. This is mostly there for gif generation (see Visuals).

The :code:`max_iters` tag is the maximum number of descent iterations allowed.

The :code:`es_crit` and :code:`dd_crit` tags are termination conditions based on the ergodic score and directional derivative, respectively. A lot of research suggests using a directional derivative criterion.

An example using some of these tags is shown below:
::
    pto_trajectory(em, tm, max_iters=1000, dd_crit=1e-4, verbose=false)
