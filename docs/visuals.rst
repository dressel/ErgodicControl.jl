=========================
Visuals
=========================

Plotting
===========
To plot the spatial distribution of an ErgodicManager, just use plot.
::

    em = ErgodicManager("single gaussian")
    plot(em)

Optional arguments customize the distribution's appearance. The code :code:`alpha` argument controls opaqueness; when :code:`alpha=1`, the distribution is fully opaque and when :code:`alpha=0`, the distribution is fully transparent. The default value is 1.

The :code:`cmap` argument controls the colormap used. The default value is :code:`"Greys"`, which uses black for the densest areas and white for the least dense areas. An alternative is the reverse, :code:`"Greys_r"`, which uses white for the densest areas and black elsewhere.

An example using these optional arguments is shown below:
::
    plot(em, cmap="Greys_r", alpha=0.8)

Likewise, you can plot a trajectory too:
::

    plot(em, xd)

In addition to the optional arguments above, there are some additional optional arguments for the trajectory:


GIFs
===========
You can create GIFs of the trajectory generation process. This allows you to see the trajectory at the end of each iteration as it is iteratively improved into the final version. To make such a GIF, you must enable logging when generating the trajectory. You can then call the :code:`gif` function on the :code:`ErgodicManager` and :code:`TrajectoryManager` used. The optional keyword argument :code:`fps` controls the frames per second. I've found 17 to be enjoyable but you can change it if you want.
::
    
    xd, ud = pto_trajectory(em, tm, max_iters=100, logging=true)
    gif(em, tm, fps=17)

You can also make a GIF of the final trajectory, which shows the agent moving along the trajectory. In most cases, this is not particularly interesting, as a simple line can capture the trajectory. However, it is more intersting if the spatial distribution evolves with time (see examples). The function call requires only the :code:`ErgodicManager` and the trajectory; the optional :code:`logging` keyword is not needed while generating the trajectory. You can also set the optional keyword argument :code:`fps` here. It is defaulted to 5.
::
    
    gif(em, xd)
