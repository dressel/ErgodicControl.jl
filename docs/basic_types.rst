=========================
Basic Types
=========================


Trajectory
=========================


Domain
=========================
The :code:`Domain` type has the following fields:
::

	# user provided
	mins::Vector{Float64}			# minimum value per dimension
	maxes::Vector{Float64}			# maximum value per dimension
	cells::Vector{Int}				# number of cells per dimension

	# calculated and used internally
	lengths::Vector{Float64}		# maxes - mins
	cell_lengths::Vector{Float64}	# length of cell in each dimension
	num_dims::Int					# number of dimensions
	cell_size::Float64				# size (area, volume, etc) of a cell

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

You should note that for :math:`\mathbb{R}^n`, the minimums need to be zero for the cosine Fourier ergodic metric to compute correctly.

To generate domains in SE(2), just ensure there are three dimensions and that the last one covers the entire angular space. The minimums need not be zero for an SE(2) domain.
::
    
    d = Domain([-1,-1,-pi], [1,1,pi], 50)

When computing ergodic trajectories over SE(2), it is recommended that you use more than 50 (or fewer) cells per dimension because SE(2) trajectories use Julia's :math:`besselj` function, which makes computation slow.


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

To learn more about plotting, check out the `Visuals` section of this readme.
