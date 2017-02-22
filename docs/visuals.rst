=========================
Visuals
=========================

Plotting
===========
To plot the spatial distribution of an ErgodicManager, just use plot.
::

    em = ErgodicManager("single gaussian")
    plot(em)

Likewise, you can plot a trajectory too:
::

    plot(em, xd)


GIFs
===========
You can create gifs of ergodic trajectory generation. Below is an example of how to generate gifs and change the frame rate.
::
    
    xd, ud = new_trajectory(em, tm, max_iters=100, logging=true)
    gif(em, tm, fps=17)
