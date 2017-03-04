######################################################################
# domain.jl
#
# experimental, not really in use yet
######################################################################
export Domain
type Domain
	mins::Vector{Float64}
	maxes::Vector{Float64}
	cells::Vector{Int}
	lengths::Vector{Float64}
	cell_lengths::Vector{Float64}
	num_dims::Int
	cell_size::Float64

	function Domain(mins::VF, maxes::VF, cells::Vector{Int})
		println("Making domain!!!")
		# TODO: check that all the dimensions match
		num_dims = length(mins)
		@assert num_dims == length(mins)
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
