######################################################################
# trajectory_generation.jl
######################################################################

#include("lqr.jl")
#include("lq.jl")
include("gradients.jl")
include("scoring.jl")
include("printing.jl")
include("new_trajectory.jl")
include("clerc_trajectory.jl")
include("cerc_trajectory.jl")
include("max_trajectory.jl")
include("kmeans_trajectory.jl")

# returns xdn and udn, the feasible projected trajectory
function project(em::ErgodicManager, tm::TrajectoryManager, K::VMF, xd::VVF, ud::VVF, zd::VVF, vd::VVF, step_size::Float64)
	xdn = [xd[1]]
	udn = Array(Vector{Float64}, 0)

	# perform descent
	alpha = VVF(tm.N + 1)
	for n = 0:tm.N
		alpha[n+1] = xd[n+1] + step_size * zd[n+1]
	end

	# perform the projection
	for n = 1:tm.N
		push!(udn, ud[n] + step_size*vd[n] + K[n]*(alpha[n] - xdn[n]))
		push!(xdn, forward_euler(tm, xdn[n], udn[n]) )
	end
	return xdn, udn
end


# modifies (xd,ud) by moving step_size in direction (zd,vd)
function descend!(xd::VV_F, ud::VV_F, zd::Matrix{Float64}, vd::Matrix{Float64}, step_size::Float64, N::Int)
	num_u = length(ud[1])
	num_x = length(xd[1])
	for i = 0:(N-1)
		for j = 1:num_u
			ud[i+1][j] += step_size*vd[j,i+1]
		end
		for j = 1:num_x
			xd[i+1][j] += step_size*zd[j,i+1]
		end
	end
	for j = 1:num_x
		xd[N+1][j] += step_size*zd[j,N+1]
	end
end

function descend!(xd::VV_F, ud::VV_F, zd::VV_F, vd::VV_F, step_size::Float64, N::Int)
	num_u = length(ud[1])
	num_x = length(xd[1])
	for i = 0:(N-1)
		for j = 1:num_u
			ud[i+1][j] += step_size*vd[i+1][j]
		end
		for j = 1:num_x
			xd[i+1][j] += step_size*zd[i+1][j]
		end
	end
	for j = 1:num_x
		xd[N+1][j] += step_size*zd[N+1][j]
	end
end


function check_convergence(es::Float64, es_crit::Float64, i::Int, max_iters::Int, dd::Float64, dd_crit::Float64, verbose::Bool, es_count::Int)
	not_finished = true
	if es < es_crit
		not_finished = false
		if verbose
			println("reached ergodic criterion...")
		end
	end
	if i > max_iters
		not_finished = false
		if verbose
			println("max iterations reached...")
		end
	end
	if abs(dd) < dd_crit
		not_finished = false
		if verbose
			println("reached directional derivative criterion...")
		end
	end
	if es_count > 50
		not_finished = false
		if verbose
			println("We've been stuck for 50 iterations...")
		end
	end
	return not_finished
end

# called if logging, not meant for general use
function save(outfile::IOStream, xd::VVF)
	n = length(xd[1])
	for xi in xd
		for i = 1:(n-1)
			wi = xi[i]
			write(outfile,"$(xi[i]),")
		end
		write(outfile,"$(xi[n])\n")
	end
end
