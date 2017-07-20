######################################################################
# domain.jl
#
# TODO: allow creation from tuples as well 
######################################################################
export Domain

type Domain
	# user provided
	mins::Vector{Float64}			# minimum value per dimension
	maxes::Vector{Float64}			# maximum value per dimension
	cells::Vector{Int}				# number of cells per dimension

	# calculated and used internally
	lengths::Vector{Float64}		# maxes - mins
	cell_lengths::Vector{Float64}	# length of cell in each dimension
	num_dims::Int					# number of dimensions
	cell_size::Float64				# size (area, volume, etc) of a cell

	function Domain(mins, maxes, cells::Vector{Int})
		num_dims = length(mins)
		@assert num_dims == length(maxes)
		@assert num_dims == length(cells)
		@assert num_dims > 0

		d = new()

		# these things stay the same
		d.mins = deepcopy(mins)
		d.maxes = deepcopy(maxes)
		d.cells = deepcopy(cells)

		# these can be easily computed
		d.num_dims = num_dims
		d.lengths = zeros(num_dims)
		d.cell_lengths = zeros(num_dims)
		d.cell_size = 1.0
		for i = 1:num_dims
			d.lengths[i] = d.maxes[i] - d.mins[i]
			d.cell_lengths[i] = d.lengths[i] / d.cells[i]
			d.cell_size *= d.cell_lengths[i]
		end
		return d
	end

	# can just pass an int for num_cells too
	function Domain(mins, maxes, num_cells::Int)
		n = length(mins)
		return Domain(mins, maxes, num_cells*ones(Int, n))
	end

	# Allow user to only give maxes
	function Domain(maxes, bins::Vector{Int}) 
		return Domain(zeros(length(maxes)), maxes, bins)
	end

	function Domain(maxes, num_cells::Int) 
		return Domain(zeros(length(maxes)), maxes, num_cells)
	end

	# allow user to only specify discretization (unit cell is assumed)
	Domain(num_cells::Int) = Domain([0,0], [1,1], [num_cells,num_cells])

end


x_min(domain::Domain) = domain.mins[1]
y_min(domain::Domain) = domain.mins[2]
z_min(domain::Domain) = domain.mins[3]

x_max(domain::Domain) = domain.maxes[1]
y_max(domain::Domain) = domain.maxes[2]
z_max(domain::Domain) = domain.maxes[3]

x_size(domain::Domain) = domain.cell_lengths[1]
y_size(domain::Domain) = domain.cell_lengths[2]
z_size(domain::Domain) = domain.cell_lengths[3]

x_cells(domain::Domain) = domain.cells[1]
y_cells(domain::Domain) = domain.cells[2]
z_cells(domain::Domain) = domain.cells[3]
