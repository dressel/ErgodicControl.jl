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

The :code:`cmap` argument controls the colormap used. The default value is :code:`Greys`, which uses black for the densest areas and white for the least dense areas. An alternative is the reverse, "Greys_r", which uses white for the densest areas and black elsewhere.

An example using these optional arguments is shown below:
::
    plot(em, cmap="Greys_r", alpha=0.8)

Likewise, you can plot a trajectory too:
::

    plot(em, xd)

In addition to the optional arguments above, there are some additional optional arguments for the trajectory:


GIFs
===========
You can create gifs of ergodic trajectory generation. Below is an example of how to generate gifs and change the frame rate.
::
    
    xd, ud = new_trajectory(em, tm, max_iters=100, logging=true)
    gif(em, tm, fps=17)
