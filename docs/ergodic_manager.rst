=========================
Ergodic Manager
=========================

The ErgodicManager type keeps track of the spatial distribution.

Fields
=========
The ErgodicManager type has the following fields:
::

    K::Int
    bins::Int
    L::Float64
    cell_size::Float64
    hk::Matrix{Float64}
    phi::Matrix{Float64}        # spatial distribution
    phik::Matrix{Float64}       # spatial distribution Fourier coefficients
    Lambdak::Matrix{Float64}
    kpixl::Matrix{Float64}

Construction
=============
The constructor for the ErgodicManager types is as follows:
::

    ErgodicManager(L::Float64, K::Int, bins::Int)


Updating Spatial Distribution
==============================
The function `phik!` can be used to update the spatial distribution `phi` and Fourier coefficients `phik`:
::

    phik!(em::ErgodicManager, d::Matrix{Float64})


Reconstrucing Spatial Distribution
====================================
Sometimes we want to reconstruct a spatial distribution from the Fourier coefficients, to see how well the Fourier coefficients capture the distribution.
::

    phi = reconstruct(em::ErgodicManager)

If you have your own set of coefficients `ck`, you can use that instead of `em.phik`:
::

    phi = reconstruct(em::ErgodicManager, ck::Matrix{Float64})



Example Managers
=================
I provide a number of pre-made ergodic managers that correspond to frequently used example domains/distributions.
::

    ErgodicManager(example_name::String; K::Int=5, bins::Int=100)

Currently, there are two valid values for example_name: "single gaussian" and "double gaussian". For exmaple, you could run:
::

    em = ErgodicManager("single gaussian", K=5, bin=100)

