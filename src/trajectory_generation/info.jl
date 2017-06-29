######################################################################
# info.jl
#
# These functions used to be in scoring but they are unrelated
######################################################################

export collect_info, optimal_info

"""
`collect_info(em::ErgodicManager, traj::VVF, d_rate::Float64)`

`collect_info(em::,rgodicManager, traj::VVF)`

modifies em.phi according to some submodular.
we don't use the last point in the trajectory

decreases at rate `D/T`
If you spend `h` time there, it is equivalent to `h*D/(h*N) = D/N`

Returns total info picked up (a scalar value).
"""
function collect_info(em::ErgodicManager, traj::VVF; steps=0)
	N = length(traj) - 1
	D = sum(em.phi)
	d_rate = D/(N+1)

	total_info = 0.0
	if steps != 0
		N = steps
	end

	for n in 0:N
		# determine the cell that this point is in
		xi,yi = find_cell(em.domain, traj[n+1])

		# if there is enough info, grab it
		info_value = min(em.phi[xi,yi], d_rate)
		em.phi[xi,yi] -= info_value
		total_info += info_value
	end
	return total_info
end



function find_cell(em::ErgodicManager, x::Vector{Float64})
	return find_cell(em.bins, em.cell_size, x)
end
# should be the proper one
function find_cell(d::Domain, x::Vector{Float64})
	x1 = round(Int, x[1] / d.cell_lengths[1], RoundDown) + 1
	x2 = round(Int, x[2] / d.cell_lengths[2], RoundDown) + 1
	if x1 > d.cells[1]; x1 -= 1; end
	if x2 > d.cells[2]; x2 -= 1; end
	if x1 < 1; x1 += 1; end
	if x2 < 1; x2 += 1; end
	return x1, x2
end


"""
`optimal_info(em, N)`

Returns optimal information gathered
"""
function optimal_info(em::ErgodicManager, N::Int)
	d_rate = sum(em.phi) / (N+1)
	total_info = 0.0

	# for each point in the trajectory...
	for n = 1:(N+1)
		best_ind = indmax(em.phi)
		info_value = min(em.phi[best_ind], d_rate)
		em.phi[best_ind] -= info_value
		total_info += info_value
	end

	return total_info
end


# create trajectory to do the above
# TODO: this should maybe go somewhere else
export optimal_traj
function optimal_traj(em::ErgodicManager, tm::TrajectoryManager)
	d_rate = sum(em.phi)/tm.N
	num_cells = em.bins*em.bins
	total_info = 0.0
	xd = Array(Vector{Float64}, tm.N+1)
	xd[1] = deepcopy(tm.x0)
	temp_phi = deepcopy(em.phi)
	size_tuple = (em.bins, em.bins)
	for n = 1:tm.N
		bi = indmax(temp_phi)
		xi, yi = ind2sub(size_tuple, bi)
		xd[n+1] = [(xi-0.5)*em.cell_size, (yi-0.5)*em.cell_size]
		temp_phi[bi] -= min(temp_phi[bi], d_rate)
	end
	ud = compute_controls(xd, tm.h)
	return xd,ud
end
