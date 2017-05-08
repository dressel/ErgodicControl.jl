######################################################################
# kmeans_trajectory.jl
######################################################################

function kmeans_trajectory(em::ErgodicManager, tm::TrajectoryManager)
	x0 = (tm.x0[1], tm.x0[2])
	return kmeans_trajectory(em, x0, tm.h, tm.N)
end

function kmeans_trajectory(em::ErgodicManager, x0::T2F, h::Float64, N::Int)
	# points is a vv_f
	points, ud0 = random_trajectory(N, h, x0)

	# show what the score was before we did the optimization
	ts = ergodic_score(em, points)
	println("ts (before) = ", round(ts,3))

	# create assignments and means matrix
	assignments = Array(Int, em.bins, em.bins)
	means = deepcopy(points)

	# TODO: select a termination condition...
	for i = 1:200
		assign_step(em, assignments, means)
		mean_step(em, assignments, means)
	end
	return means, ud0

end

function assign_step(em::ErgodicManager, assignments::Matrix{Int}, means::VVF)
	N = length(means)
	for xi = 1:em.bins
		for yi = 1:em.bins
			best_n = 1
			best_d = Inf
			for n = 1:N
				d = euc_dist_2(em, xi, yi, means[n])
				if d < best_d
					best_d = d
					best_n = n
				end
			end
			assignments[xi,yi] = best_n
		end
	end
end

# compute the mean of each cluster
function mean_step(em::ErgodicManager, assignments::Matrix{Int}, means::VVF)
	N = length(means)
	# let's clear the means vector
	for i = 1:N
		means[i][1] = 0.
		means[i][2] = 0.
	end

	# TODO: must we remake this every time
	mean_counts = zeros(N)

	for xi = 1:em.bins
		for yi = 1:em.bins
			i = assignments[xi,yi]
			phi_q = em.phi[xi,yi]
			means[i][1] += phi_q*(xi-0.5)*em.cell_size
			means[i][2] += phi_q*(yi-0.5)*em.cell_size
			mean_counts[i] += phi_q
		end
	end

	# normalize the means
	for i = 1:N
		mc = mean_counts[i]
		means[i][1] /= mc
		means[i][2] /= mc
	end
end

# returns the squared euclidean distance to some stuff
# TODO: give this a respectable name
function euc_dist_2(em::ErgodicManager, xi::Int, yi::Int,p::Vector{Float64})
	return euc_dist_2(em, xi, yi, p[1], p[2])
end
function euc_dist_2(em::ErgodicManager, xi::Int, yi::Int, xf::Float64, yf::Float64)
	xp = em.cell_size * (xi - 0.5)
	yp = em.cell_size * (yi - 0.5)
	dx = xf - xp
	dy = yf - yp

	return dx*dx + dy*dy
end

function random_trajectory(N::Int, h::Float64, x0::T2F)
	xd = [[x0[1], x0[2]] for i = 1:N+1]
	points = Array(Vector{Float64}, N)
	for i = 1:N
		points[i] = [rand(), rand()]
	end

	tsp_nn!(xd, points)

	ud = compute_controls(xd, h)

	return xd, ud
end
