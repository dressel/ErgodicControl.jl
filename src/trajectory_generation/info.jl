######################################################################
# info.jl
#
# These functions used to be in scoring but they are unrelated
######################################################################

"""
`collect_info(em::ErgodicManager, traj::VVF, d_rate::Float64)`

`collect_info(em::ErgodicManager, traj::VVF)`

modifies em.phi according to some submodular.
we don't use the last point in the trajectory

decreases at rate `D/T`
If you spend `h` time there, it is equivalent to `h*D/(h*N) = D/N`

Returns total info picked up (a scalar value).
"""
function collect_info(em::ErgodicManager, traj::VVF; steps=0, right::Bool=true)
	N = length(traj) - 1
	D = sum(em.phi)
	d_rate = D/N
	collect_info(em.phi, em.cell_size, traj,d_rate,steps=steps, right=right)
end
function collect_info(em::ErgodicManager, traj::VVF, d_rate::Float64; steps=0, right::Bool=true)
	collect_info(em.phi, em.cell_size, traj,d_rate,steps=steps, right=right)
end

function collect_info(phi::Matrix{Float64}, cell_size::Float64, traj::VVF; steps=0, right::Bool=true)
	N = length(traj) - 1
	D = sum(phi)
	d_rate = D/N
	collect_info(phi, cell_size, traj, d_rate, steps=steps, right=right)
end
function collect_info(phi::Matrix{Float64}, cell_size::Float64, traj::VVF, d_rate::Float64; steps=0, right::Bool=true)
	bins, rar = size(phi)
	N = length(traj) - 1
	total_info = 0.0
	if steps != 0
		N = steps
	end
	N_range = 0:(N-1)
	if right
		N_range += 1
	end
	for n in N_range
		xi,yi = find_cell(bins, cell_size, traj[n+1])

		# if there is enough info, grab that shit yo
		info_value = min(phi[xi,yi], d_rate)
		phi[xi,yi] -= info_value
		total_info += info_value
	end
	return total_info
end
export collect_info


function find_cell(em::ErgodicManager, x::Vector{Float64})
	return find_cell(em.bins, em.cell_size, x)
end
function find_cell(bins::Int, cell_size::Float64, x::Vector{Float64})
	x1 = round(Int, x[1] / cell_size, RoundDown) + 1
	x2 = round(Int, x[2] / cell_size, RoundDown) + 1
	if x1 > bins; x1 -= 1; end
	if x2 > bins; x2 -= 1; end
	if x1 < 1; x1 += 1; end
	if x2 < 1; x2 += 1; end
	return x1, x2
end

# 
function optimal_info(em::ErgodicManager, N::Int)
	d_rate = sum(em.phi)/N
	num_cells = em.bins*em.bins
	total_info = 0.0
	for n = 1:N
		best_i = 0
		for i = 1:num_cells
			if em.phi[i] > d_rate
				best_i = i
				total_info += d_rate
				em.phi[best_i] -= d_rate
				break
			end
		end
		if best_i == 0  # we didn't find a good enough cell,
			# loop over and find max
			best_i = indmax(em.phi)
			total_info += em.phi[best_i]
			em.phi[best_i] = 0.0
		end
	end
	return total_info
end
export optimal_info

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
