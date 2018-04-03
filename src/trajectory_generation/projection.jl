######################################################################
# projection.jl
#
# Handles projections and descents
######################################################################

# returns xdn and udn, the feasible projected trajectory
function project(em::ErgodicManager, tm::TrajectoryManager, K::VMF, xd::VVF, ud::VVF, zd::VVF, vd::VVF, step_size::Float64)
	xdn = [xd[1]]
	udn = Array{VF}(0)

	# perform descent
	alpha = VVF(tm.N + 1)
	for n = 0:tm.N
		alpha[n+1] = xd[n+1] + step_size * zd[n+1]
	end

	# perform the projection
	for n = 1:tm.N
		push!(udn, ud[n] + step_size*vd[n] + K[n]*(alpha[n] - xdn[n]))
		push!(xdn, integrate(tm, xdn[n], udn[n]) )
	end
	return xdn, udn
end

# A projection for LTI systems
function project2(em::ErgodicManager, tm::TrajectoryManager, K::VMF, xd::VVF, ud::VVF, zd::VVF, vd::VVF, step_size::Float64)
	xdn = [xd[1]]
    udn = Array{VF}(0)

	xdn = VVF(0)
	udn = VVF(0)

	# perform the projection
	# Shouldn't need to even integrate...
	for n = 1:tm.N
		push!(udn, ud[n] + step_size*vd[n])
		push!(xdn, xd[n] + step_size*zd[n])
		#push!(xdn, integrate(tm, xdn[n], udn[n]) )
		#push!(xdn, symplectic_euler(tm, xdn[n], udn[n]) )
	end
	push!(xdn, xd[tm.N+1] + step_size*zd[tm.N+1])
	return xdn, udn
end


# No projection and performed in place
function descend!(xd::VVF, ud::VVF, zd::VVF, vd::VVF, step_size::Float64)
	N = length(ud)
	for ni = 1:N
		ud[ni] += step_size * vd[ni]
		xd[ni] += step_size * zd[ni]
	end
	xd[N+1] += step_size * zd[N+1]
end
