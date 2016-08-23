######################################################################
# clerc_trajectory.jl
# CLErC (Constrained Linear Ergodic Control)
# Assumes linear dynamics, but you can constrain norm of effort
######################################################################

using Convex, SCS

function clerc_trajectory(em::ErgodicManager, tm::TrajectoryManager; verbose::Bool=true, logging::Bool=false, max_iters::Int=30)
	xd0, ud0 = initialize(tm.initializer, em, tm)
	clerc_trajectory(em, tm, xd0, ud0; verbose=verbose, logging=logging, max_iters=max_iters)
end

function clerc_trajectory(em::ErgodicManager, tm::TrajectoryManager, xd0::VV_F, ud0::VV_F; verbose::Bool=true, logging::Bool=false, max_iters::Int=30, es_crit::Float64=0.003)

	# let's not overwrite the initial trajectories
	xd = deepcopy(xd0)
	ud = deepcopy(ud0)
	N = tm.N

	# matrices for gradients
	ad = zeros(tm.n, N)
	bd = zeros(tm.m, N)

	# prepare variables for convex optimization
	z = Variable(tm.n, N+1)
	v = Variable(tm.m, N)
	c = Array(Constraint, 0)
	push!(c, z[:,1] == 0)
	for n = 0:(N-1)
		push!(c, z[:,n+1+1] == tm.A*z[:,n+1] + tm.B*v[:,n+1])
	end

	# prepare for logging if need be 
	if logging
		outfile = open("temp.csv", "a")
		save(outfile, xd)
	end

	if verbose; print_header(); end
	i = 1
	not_finished = true
	es = 0.; cs = 0.; ts = 0.; dd = 0.; step_size = 0.
	while not_finished

		# Find gradients, descent step, and step size. Then descend!
		gradients!(ad, bd, em, tm, xd, ud)
		zd, vd = convex_descent(ad, bd, N, z, v, c, tm.R[1,1])
		step_size = .15 / sqrt(i)
		#step_size = .01 / sqrt(i)
		descend!(xd, ud, zd, vd, step_size, N)

		# compute statistics and report
		es, cs, ts = all_scores(em, tm, xd, ud)
		dd = directional_derivative(ad, bd, zd, vd)
		#sdd = scaled_dd(ad, bd, zd, vd)
		if verbose; step_report(i, es, cs, ts, dd, step_size); end
		if logging; save(outfile, xd); end

		# check convergence
		i += 1
		not_finished = check_convergence(es,es_crit,i,max_iters,verbose)
	end

	# now that we are done, print a special finished report
	if verbose
		print_header()
		if verbose; step_report(i-1, es, cs, ts, dd, step_size); end
	end

	if logging; close(outfile); end

	return xd, ud
end

function check_convergence(es::Float64, es_crit::Float64, i::Int, max_iters::Int, verbose::Bool)
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
	return not_finished
end

function step_report(i::Int, es::Float64, cs::Float64, ts::Float64, dd::Float64, step_size::Float64)
	@printf " %-7i" i
	@printf " %-14.7f" es
	@printf " %-14.7f" cs
	@printf " %-12.7f" ts
	@printf " %-12.7f" dd
	@printf " %-7.5f" step_size
	println()
end

function print_header()
	print_dashes()
	println(" iter  |ergodic score |control score |total score |direc deriv |step size")
	print_dashes()
end
function print_dashes()
	println("--------------------------------------------------------------------------")
end

function gradients!(ad::Matrix{Float64}, bd::Matrix{Float64}, em::ErgodicManager, tm::TrajectoryManager, xd::VV_F, ud::VV_F)
	for n = 0:(tm.N-1)
		an_x, an_y = compute_ans(em, xd, tm, n)
		ad[1,n+1] = an_x
		ad[2,n+1] = an_y
		# assume that only the first two state variables matter
		# in fact, this isn't necessary if it's been initalized to zeros...
		for i = 3:tm.n
			ad[i,n+1] = 0.0
		end

		#bd[n+1] = tm.h * tm.R * ud[n+1]
		# TODO: allow for different tm.m
		bd[1,n+1] = tm.h * (tm.R[1,1]*ud[n+1][1] + tm.R[1,2]*ud[n+1][2])
		bd[2,n+1] = tm.h * (tm.R[2,1]*ud[n+1][1] + tm.R[2,2]*ud[n+1][2])
	end
end

function compute_ans(em::ErgodicManager, xd::VV_F, tm::TrajectoryManager, n::Int)
	xnx = xd[n+1][1]
	xny = xd[n+1][2]
	L = em.L

	an_x = 0.0
	an_y = 0.0
	 
	for k1 = 0:em.K
		for k2 = 0:em.K
			hk = em.hk[k1+1,k2+1]

			dFk_dxn1 = -k1*pi*sin(k1*pi*xnx/L)*cos(k2*pi*xny/L) / (hk*L)
			dFk_dxn2 = -k2*pi*cos(k1*pi*xnx/L)*sin(k2*pi*xny/L) / (hk*L)

			fk = 0.0
			for i = 0:(tm.N-1)
				x = xd[i+1]
				fk += cos(k1*pi*x[1]/L) * cos(k2*pi*x[2]/L) / hk
			end
			c = em.Lambdak[k1+1,k2+1] * (tm.h*fk/tm.T - em.phik[k1+1,k2+1])
			an_x += c*dFk_dxn1
			an_y += c*dFk_dxn2
		end
	end
	an_x *= 2.0*tm.h
	an_y *= 2.0*tm.h
	return an_x, an_y
end

# descent direction using convex optimization
function convex_descent(a::Matrix{Float64}, b::Matrix{Float64}, N::Int, z::Variable, v::Variable, c::Vector{Constraint}, r::Float64)

	# create the problem
	problem = minimize(vecdot(a,z[:,1:N]) + vecdot(b,v) + sumsquares(z) + r*sumsquares(v), c)

	# solve the problem
	solve!(problem, SCSSolver(verbose=0,max_iters=100000))

	# put z and v into vec of vec format
	#vd = Array(Vector{Float64}, N+1)
	#zd = Array(Vector{Float64}, N+1)
	#for i = 1:N
	#	zd[i] = [z.value[1,i], z.value[2,i]]
	#	vd[i] = [v.value[1,i], v.value[2,i]]
	#end
	## just repeat last control, it technically doesn't matter
	#zd[N+1] = [z.value[1,N+1], z.value[2,N+1]]
	#vd[N+1] = [v.value[1,N], v.value[2,N]]

	#return zd, vd
	return z.value, v.value
end


# modifies (xd,ud) by moving step_size in direction (zd,vd)
function descend!(xd::VV_F, ud::VV_F, zd::Matrix{Float64}, vd::Matrix{Float64}, step_size::Float64, N::Int)
	for i = 0:(N-1)
		ud[i+1][1] += step_size*vd[1,i+1]
		ud[i+1][2] += step_size*vd[2,i+1]
		xd[i+1][1] += step_size*zd[1,i+1]
		xd[i+1][2] += step_size*zd[2,i+1]
	end
	ud[N+1] = ud[N]
end

function save(outfile::IOStream, xd::VV_F)
	n = length(xd[1])
	for xi in xd
		for i = 1:(n-1)
			wi = xi[i]
			write(outfile,"$(xi[i]),")
		end
		write(outfile,"$(xi[n])\n")
	end
end
