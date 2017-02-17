=========================
Visuals
=========================

Plotting
===========
To plot the spatial distribution of an ErgodicManager, just use plot.
::

    em = ErgodicManager("single gaussian")
    plot(em)

GIFs
===========
You can create gifs of ergodic trajectory generation.

Below is an example of how to generate gifs. This is ugly and will be updated in the future.
::
    
    xd, ud = new_trajectory(em, tm, max_iters=100, logging=true)
    trajectories = readcsv("temp.csv")
    gif(em, trajectories, 101, fps=17)
