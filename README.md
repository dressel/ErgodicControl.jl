# ErgodicControl

This package currently generates ergodic trajectories for linear systems.
Currently, the distributions that ergodic trajectories are made ergodic to must be 2-dimensional and square.

In the future, I hope to add support for non-linear systems and adding additional dimensions (at least getting to SE(2)).

## ErgodicManager
The `ErgodicManager` type handles a lot of the stuff we care about.
It stores the distribution `phi`, the coefficients `phi_k`, etc.
It can be constructed by providing the side length of a square `L`, the number of coefficients `K`, and the number of bins per side length used in discretization:
```
ErgodicManager(L::Float64, K::Int, bins::Int)
```
There are also two preloaded examples that can be loaded by giving an input string.
Valid values of `example_name` are "single gaussian" and "double gaussian".
```
ErgodicManager(example_name::ASCIIString; K::Int=5, bins::Int=100)
```

## TrajectoryManager
Before generating trajectories, you need to create a `TrajectoryManager`, which stores a lot of information used during trajectory generation.

## Generating Trajectories


## Plotting and GIF Generation
Plotting is done with PyPlot.jl.

[![Build Status](https://travis-ci.org/dressel/ErgodicControl.jl.svg?branch=master)](https://travis-ci.org/dressel/ErgodicControl.jl)
