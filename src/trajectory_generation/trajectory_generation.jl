######################################################################
# trajectory_generation.jl
######################################################################

#include("lqr.jl")
#include("lq.jl")
include("scoring.jl")
include("new_trajectory.jl")
include("clerc_trajectory.jl")
include("cerc_trajectory.jl")
include("max_trajectory.jl")
include("kmeans_trajectory.jl")

# TODO: do dynamics correctly
# returns xdn and udn, the feasible projected trajectory
function project(em::ErgodicManager, tm::TrajectoryManager, K::Vector{MF}, xd::VV_F, ud::VV_F, zd::VV_F, vd::VV_F, step_size::Float64)
	xdn = [xd[1]]
	udn = Array(Vector{Float64}, 0)
	for n = 1:tm.N
		push!(udn, ud[n] + step_size*vd[n] + step_size*K[n]*zd[n])
		push!(xdn, tm.dynamics.A*xdn[n] + tm.dynamics.B*udn[n])
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
