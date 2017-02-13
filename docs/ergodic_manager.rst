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
    phi::Matrix{Float64}
    phik::Matrix{Float64}
    Lambdak::Matrix{Float64}
    kpixl::Matrix{Float64}

Construction
=============
The constructor for the ErgodicManager types is as follows:
::
ErgodicManager(L::Float64, K::Int, bins::Int)


Example Managers
=================
I provide a number of pre-made ergodic managers that correspond to frequently used example domains/distributions.
::

    ErgodicManager(example_name::String; K::Int=5, bins::Int=100)

Currently, there are two valid values for example_name: "single gaussian" and "double gaussian". For exmaple, you could run:
::
em = ErgodicManager("single gaussian, K=5, bin=100)

