######################################################################
# scoring.jl
# Computes ergodic, total, and control scores
######################################################################

"""
`decompose(em, traj::VV_F, start_idx=1)`

Decomposes a set of positions into a set of `ck` Fourier coefficients.

If `start_idx` is 1, then the right Riemann sum is used.
If `start_idx` is 0, then the left Riemann sum is used.
"""
#function decompose(em::ErgodicManager, traj::VV_F, start_idx::Int=0)
#	traj2 = [(traj[i][1], traj[i][2]) for i = 1:length(traj)]
#	return decompose(em, traj2, start_idx)
#end
#function decompose(em::ErgodicManager, traj::V_T2F)
#	K = em.K
#	N = length(traj)-1
#	ck = zeros(K+1, K+1)
#	for k1 = 0:K
#		kpiL1 = k1 * pi / em.L
#		for k2 = 0:K
#			kpiL2 = k2 * pi / em.L
#			hk = em.hk[k1+1, k2+1]
#			fk_sum = 0.0
#			# now loop over time
#			for n = 0:N-1
#				xn = traj[n+1]
#				fk_sum += cos(kpiL1 * xn[1])  * cos(kpiL2 * xn[2])
#			end
#			ck[k1+1, k2+1] = fk_sum / (hk * N)
#		end
#	end
#	return ck
#end


#function decompose(em::ErgodicManager, traj::V_T2F, start_idx::Int=0)
function decompose(em::ErgodicManager, traj::VV_F, start_idx::Int=1)
	K = em.K
	N = length(traj)-1
	ck = zeros(K+1, K+1)
	for k1 = 0:K
		kpiL1 = k1 * pi / em.L
		for k2 = 0:K
			kpiL2 = k2 * pi / em.L
			hk = em.hk[k1+1, k2+1]
			fk_sum = 0.0
			# now loop over time
			for n = 0:N-1
				xn = traj[n + start_idx + 1]
				fk_sum += cos(kpiL1 * xn[1])  * cos(kpiL2 * xn[2])
			end
			ck[k1+1, k2+1] = fk_sum / (hk * N)
		end
	end
	return ck
end

function reconstruct(em::ErgodicManager, xd::VV_F, start_idx::Int=1)
	ck = decompose(em, xd, start_idx)
	return reconstruct(em, ck)
end


"""
`ergodic_score(em, traj::V_T2F)`

First breaks down the trajectory into components ck.
"""
function ergodic_score(em::ErgodicManager, traj::VV_F, start_idx::Int=1)
	ck = decompose(em, traj, start_idx)
	return ergodic_score(em, ck)
end
function ergodic_score(em::ErgodicManager, ck::Matrix{Float64})
	val = 0.0
	for k1 = 0:em.K
		for k2 = 0:em.K
			d = em.phik[k1+1,k2+1] - ck[k1+1,k2+1]
			val += em.Lambdak[k1+1,k2+1] * d * d
		end
	end
	return val
end

"""
`control_score(ud::VV_F, R, h)`

Assumes only non-zero elements of `R` are corners.

`control_score(ud::VV_F, N::Int)`

Computes sum_i=1^N dot(u,u). Does not divide by 2 or do anything else.
"""
function control_score(ud::VV_F, R::Matrix{Float64}, h::Float64)
	cs = 0.0
	for ui in ud
		cs += R[1,1] * ui[1] * ui[1]
		cs += R[2,2] * ui[2] * ui[2]
	end
	return 0.5 * h * cs
end
control_score(ud::VV_F) = control_score(ud, eye(2), 1.0)
function control_score(ud::VV_F, N::Int)
	cs = 0.0
	for n = 1:N
		cs += ud[n][1] * ud[n][1]
		cs += ud[n][2] * ud[n][2]
	end
	return cs
end

function control_effort(ud::VV_F, N::Int)
	cs = 0.0
	for n = 1:N
		udx = ud[n][1]
		udy = ud[n][1]
		cs += sqrt(udx*udx + udy*udy)
	end
	return cs
end


"""
`total_score(em, xd::VV_F, ud::VV_F, T::Float64)`

Computes the total score `q*ergodic_score + sum_n h/2 un'Rn un`

Currently assumes `q = 1.0` and `R = 0.01 * eye(2)`
"""
function total_score(em::ErgodicManager, tm::TrajectoryManager, xd::VV_F, ud::VV_F)
	return q * ergodic_score(em, xd) + control_score(ud, tm.R, tm.h)
end
# TODO: actually get q and R from the correct place 
function total_score(em::ErgodicManager, xd::VV_F, ud::VV_F, T::Float64)
	q = 1.0
	R = 0.01 * eye(2)
	N = length(xd) - 1
	h = T/N
	return q * ergodic_score(em, xd) + control_score(ud, R, h)
end
# TODO: let's not make this so shitty...
function total_score(em::ErgodicManager, xd::VV_F, ud::VV_F, zd::VV_F, vd::VV_F, alpha::Float64, T::Float64)
	xd2 = deepcopy(xd)
	ud2 = deepcopy(ud)
	for i = 1:length(xd2)
		xd2[i][1] += alpha * zd[i][1]
		xd2[i][2] += alpha * zd[i][2]

		ud2[i][1] += alpha * vd[i][1]
		ud2[i][2] += alpha * vd[i][2]
	end
	return total_score(em, xd2, ud2, T)
end

"""
`all_scores(em::ErgodicManager, tm::TrajectoryManager, xd::VV_F, ud::VV_F)`
"""
#function all_scores(em::ErgodicManager, tm::TrajectoryManager, xd::VV_F, ud::VV_F, N_range::UnitRange{Int})
#	es = ergodic_score(em, xd, N_range)
#	cs = control_score(ud, tm.R, tm.h)
#	ts = tm.q*es + cs
#	return es, cs, ts
#end
function all_scores(em::ErgodicManager, tm::TrajectoryManager, xd::VV_F, ud::VV_F, start_idx::Int=1)
	es = ergodic_score(em, xd, start_idx)
	cs = control_score(ud, tm.R, tm.h)
	ts = tm.q*es + cs
	return es, cs, ts
end


"""
`collect_info(em::ErgodicManager, traj::VV_F, d_rate::Float64)`

`collect_info(em::ErgodicManager, traj::VV_F)`

modifies em.phi according to some submodular.
we don't use the last point in the trajectory

decreases at rate `D/T`
If you spend `h` time there, it is equivalent to `h*D/(h*N) = D/N`

Returns total info picked up (a scalar value).
"""
function collect_info(em::ErgodicManager, traj::VV_F; steps=0, right::Bool=true)
	N = length(traj) - 1
	D = sum(em.phi)
	d_rate = D/N
	collect_info(em.phi, em.cell_size, traj,d_rate,steps=steps, right=right)
end
function collect_info(em::ErgodicManager, traj::VV_F, d_rate::Float64; steps=0, right::Bool=true)
	collect_info(em.phi, em.cell_size, traj,d_rate,steps=steps, right=right)
end

function collect_info(phi::Matrix{Float64}, cell_size::Float64, traj::VV_F; steps=0, right::Bool=true)
	N = length(traj) - 1
	D = sum(phi)
	d_rate = D/N
	collect_info(phi, cell_size, traj, d_rate, steps=steps, right=right)
end
function collect_info(phi::Matrix{Float64}, cell_size::Float64, traj::VV_F, d_rate::Float64; steps=0, right::Bool=true)
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
	if x1 < 0; x1 += 0; end
	if x2 < 0; x2 += 0; end
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
