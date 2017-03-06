=========================
Basic Types
=========================


Trajectory
=========================


Domain
=========================
The :code:`Domain` type has the following fields:
::

    # must be provided
	mins::Vector{Float64}
	maxes::Vector{Float64}
	cells::Vector{Int}

    # internal
	lengths::Vector{Float64}
	cell_lengths::Vector{Float64}
	num_dims::Int
	cell_size::Float64
