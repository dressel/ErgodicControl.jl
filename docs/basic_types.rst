=========================
Basic Types
=========================


Type Aliases
=========================
Type aliases are used to simplify tedious code. For example, the type corresponding to a vector of vector of floats is :code:`Vector{Vector{Float64}}`. To make function arguments easier to read, this is type aliased to the simpler :code:`VVF`. Below is a list of type aliases used throughout the code and documentation.
::

    typealias MF  Matrix{Float64}
    typealias VMF Vector{MF}
    typealias VF  Vector{Float64}
    typealias VVF Vector{VF}


Trajectory
=========================
A discrete trajectory is simply a set of states. Each state is represented a a vector. For example, a point in 2D space is represented with a vector of length 2.

A trajectory is a vector of state vectors (so :code:`VVF` according to the aliases above). 
This might seem tedious, and a simpler representation might have been a 2D array, where each row represents a single state.
The reasoning behind the :code:`VVF` representation is that slicing vectors out of a 2D array slows Julia code down.
I'm not sure if it has a noticeable effect (if any), but that was the thought process.

To convert between these representations, the :code:`traj2mat` and :code:`mat2traj` functions have been included.


Domain
=========================
The :code:`Domain` type has the following fields:
::

	# user provided
	mins::Vector{Float64}           # minimum value per dimension
	maxes::Vector{Float64}          # maximum value per dimension
	cells::Vector{Int}              # number of cells per dimension

	# calculated and used internally
	lengths::Vector{Float64}        # maxes - mins
	cell_lengths::Vector{Float64}   # length of cell in each dimension
	num_dims::Int                   # number of dimensions
	cell_size::Float64              # size (area, volume, etc) of a cell

The :code:`Domain` constructor requires knowledge of the domain minimums, maximums, and number of cells per side.
::
    
    Domain(mins, maxes, cells)

If you wanted to create the unit square with each dimension discretized into 100 bins (100^2 for the whole square), you have the following options:
::

    # verbose
    d = Domain([0,0], [1,1], [100,100])

    # if one discretization level is provided, all dimensions assume it
    d = Domain([0,0], [1,1], 100)

    # if you don't provide minimums, they are assumed to be zero
    d = Domain([1,1], [100,100])
    d = Domain([1,1], 100)

    # Unit square in R^2 with 100 cells per side
    d = Domain(100)

To generate domains in SE(2), just ensure there are three dimensions and that the last one covers the entire angular space.
::
    
    d = Domain([-1,-1,-pi], [1,1,pi], 50)

When computing ergodic trajectories over SE(2), it is recommended that you use 50 (or fewer) cells per dimension because SE(2) trajectories use Julia's :code:`besselj` function, which makes computation slow.


Probability Distributions
===========================
Probability distributions over :math:`\mathbb{R}^n` are simply represented as arrays with :math:`n` dimensions.

The provided :code:`gaussian` functions make it easy to generate distributions over a domain
::
    
    d = Domain([2,1], [200,100])

    # returns array with 100 rows and 200 columns
    phi = gaussian(d, [1.5,0.5], 0.03*eye(2))

    # plot
    plot(d, phi)
    xlabel("x")
    ylabel("y")

The resulting plot is shown below. To learn more about plotting, check out the `Visuals` section of this readme.

.. image:: http://stanford.edu/~dressel/pics/single_domain.png

The :code:`gaussian` function can also do multi-Gaussians. In this case the function argument is :code:`gaussian(domain, means, Sigmas, weights)`. The :code:`means` argument is a vector of mean vectors; :code:`Sigmas` is a vector of covariance matrices; :code:`weights` is a vector of weights, one for each Gaussian. This argument is optional and the default value weights matrices equally. The :code:`gaussian` function will ensure the resulting distribution is a proper density, so you the weights don't need to add to 1.
::

    d = Domain([1,1], 100)

    means = [[.3,.7], [.7,.3]]
    Sigmas = [.025*eye(2), .025*eye(2)]
    weights = [1,2]
    phi = gaussian(d, means, Sigmas, weights)

    plot(d, phi)
    xlabel("x")
    ylabel("y")

.. image:: http://stanford.edu/~dressel/pics/multi_domain.png
