######################################################################
# scoring.jl
# Computes ergodic, total, and control scores
######################################################################

"""
`decompose(em, traj::VVF)`

Decomposes a set of positions into a set of `ck` Fourier coefficients.

"""
function decompose(em::ErgodicManagerR2, traj::VVF)
	K = em.K
	N = length(traj)-1
	#N = length(traj)
	ck = zeros(K+1, K+1)
	Lx = em.domain.lengths[1]
	Ly = em.domain.lengths[2]
	xm = x_min(em)
	ym = y_min(em)
	for k1 = 0:K
		kpiL1 = k1 * pi / Lx
		for k2 = 0:K
			kpiL2 = k2 * pi / Ly
			hk = em.hk[k1+1, k2+1]
			fk_sum = 0.0
			# now loop over time
			#for n = 0:N-1
			for n = 0:N
				xn = traj[n + 1]
				fk_sum += cos(kpiL1 * (xn[1]-xm))  * cos(kpiL2 * (xn[2]-ym))
			end
			#ck[k1+1, k2+1] = fk_sum / (hk * N)
			ck[k1+1, k2+1] = fk_sum / (hk * (N+1))
		end
	end
	return ck
end

function decompose(em::ErgodicManagerSE2, traj::VVF)
	N = length(traj)-1
	ck = zeros(Complex{Float64}, em.M+1, em.N+1, em.P+1)
	for m = 0:em.M
		for n = 0:em.N
			for p = 0:em.P
				fk_sum = 0.0im
				# now loop over time
				for i = 0:N-1
					xi = traj[i + 1]
					fk_sum += F_mnp(m,n,p,xi[1],xi[2],xi[3])
				end
				# TODO: check that this is right
				ck[m+1, n+1, p+1] = fk_sum / N
			end
		end
	end
	return ck
end

function reconstruct(em::ErgodicManager, xd::VVF)
	ck = decompose(em, xd)
	return reconstruct(em, ck)
end


"""
`ergodic_score(em, traj::VVF)`

First breaks down the trajectory into components ck.
"""
function ergodic_score(em::ErgodicManager, traj::VVF)
	ck = decompose(em, traj)
	return ergodic_score(em, ck)
end
function ergodic_score(em::ErgodicManagerR2, ck::Matrix{Float64})
	val = 0.0
	for k1 = 0:em.K
		for k2 = 0:em.K
			d = em.phik[k1+1,k2+1] - ck[k1+1,k2+1]
			val += em.Lambda[k1+1,k2+1] * d * d
		end
	end
	return val
end

function ergodic_score(em::ErgodicManagerSE2, ck::Array{Complex{Float64},3})
	val = 0.0
	for m = 0:em.M
		for n = 0:em.N
			for p = 0:em.P
				d = em.phik[m+1,n+1,p+1] - ck[m+1,n+1,p+1]
				dr = real(d)
				di = imag(d)
				val += em.Lambda[m+1,n+1,p+1] * (dr*dr + di*di)
			end
		end
	end
	return val
end

"""
`control_score(ud::VVF, R, h)`

Assumes only non-zero elements of `R` are corners (not anymore I think)

`control_score(tm::TrajectoryManager,u::VVF) = control_score(u, tm.R, tm.h)`
"""
function control_score(ud::VVF, R::Matrix{Float64}, h::Float64)
	cs = 0.0
	num_u = length(ud[1])
	for ui in ud

		for j = 1:num_u
		   cs += R[j,j] * ui[j] * ui[j]
		end

		# TODO: this is how I do it. Benchmark it
		#cs += dot(ui, R*ui)

		# old way I suppose. delete this dumb stuff
		#cs += R[1,1] * ui[1] * ui[1]
		#cs += R[2,2] * ui[2] * ui[2]
	end
	return 0.5 * h * cs
end

function control_score(tm::TrajectoryManager, ud::VVF)
	return control_score(ud, tm.R, tm.h)
end


# TODO: do I even use these functions?
control_score(ud::VVF) = control_score(ud, eye(2), 1.0)
function control_score(ud::VVF, N::Int)
	cs = 0.0
	for n = 1:N
		num_u = length(ud[n])
		for j = 1:num_u
			cs += ud[n][j] * ud[n][j]
		end
		#cs += ud[n][1] * ud[n][1]
		#cs += ud[n][2] * ud[n][2]
	end
	return cs
end

function control_effort(ud::VVF, N::Int=length(ud))
	cs = 0.0
	for n = 1:N
		udx = ud[n][1]
		udy = ud[n][1]
		cs += sqrt(udx*udx + udy*udy)
	end
	return cs
end


# if anyone is outside the domain...
#function barrier_score(em::ErgodicManager, xd::VVF, c::Float64)
#	if c == 0.0; return 0.0; end
#
#	bs = 0.0
#	for xi in xd
#		if (xi[1] > em.L) || (xi[1] < 0.0)
#			dx = (xi[1] - 0.5*em.L)
#			bs += c * dx * dx
#		end
#		if (xi[2] > em.L) || (xi[2] < 0.0)
#			dy = (xi[2] - 0.5*em.L)
#			bs += c * dy * dy
#		end
#	end
#	return bs
#end

function barrier_score(em::ErgodicManager, xd::VVF, c::Float64)
	if c == 0.0; return 0.0; end

	bs = 0.0
	xmax = x_max(em)
	ymax = y_max(em)
	xmin = x_min(em)
	ymin = y_min(em)

	for xi in xd
		if (xi[1] > xmax)
			dx = xi[1] - xmax
			bs += c * dx * dx
		elseif (xi[1] < xmin)
			dx = xi[1] - xmin
			bs += c * dx * dx
		end
		if (xi[2] > ymax)
			dy = xi[2] - ymax
			bs += c * dy * dy
		elseif (xi[2] < ymin)
			dy = xi[2] - ymin
			bs += c * dy * dy
		end
	end
	return bs
end


"""
`total_score(em::ErgodicManager, tm::TrajectoryManager, xd::VVF, ud::VVF)`

Computes the total score `q*ergodic_score + sum_n h/2 un'Rn un`
"""
function total_score(em::ErgodicManager, tm::TrajectoryManager, xd::VVF, ud::VVF)
	# TODO: add quadratic barrier score if need be
	es = tm.q * ergodic_score(em, xd)
	cs = control_score(ud, tm.R, tm.h)
	bs = barrier_score(em, xd, tm.barrier_cost)
	#bs = 0.0
	return es + cs + bs
end

"""
`all_scores(em::ErgodicManager, tm::TrajectoryManager, xd::VVF, ud::VVF)`
"""
function all_scores(em::ErgodicManager, tm::TrajectoryManager, xd::VVF, ud::VVF)
	es = ergodic_score(em, xd)
	cs = control_score(ud, tm.R, tm.h)
	ts = tm.q*es + cs
	return es, cs, ts
end


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
