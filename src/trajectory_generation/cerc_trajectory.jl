######################################################################
# cerc_trajectory.jl
# CErC (Constrained Ergodic Control)
# Does NOT assume linear dynamics
# written 2/07/2017 for Dubins stuff
######################################################################

using Convex, SCS

function cerc_trajectory(em::ErgodicManager, tm::TrajectoryManager; verbose::Bool=true, logging::Bool=false, max_iters::Int=30, es_crit::Float64=0.003,right::Bool=true)
	xd0, ud0 = initialize(tm.initializer, em, tm)
	cerc_trajectory(em, tm, xd0, ud0; verbose=verbose, logging=logging, max_iters=max_iters, es_crit=es_crit,right=right)
end

function cerc_trajectory(em::ErgodicManager, tm::TrajectoryManager, xd0::VVF, ud0::VVF; verbose::Bool=true, logging::Bool=false, max_iters::Int=30, es_crit::Float64=0.003, right::Bool=true)

	# let's not overwrite the initial trajectories
	xd = deepcopy(xd0)
	ud = deepcopy(ud0)
	N = tm.N
	start_idx = right ? 1 : 0

	# matrices for gradients
	ad = zeros(tm.dynamics.n, N)
	bd = zeros(tm.dynamics.m, N)

	# prepare variables for convex optimization
	z = Variable(tm.dynamics.n, N+1)
	v = Variable(tm.dynamics.m, N)
	c = Array(Constraint, 0)
	push!(c, z[:,1] == 0)
	#for n = 0:(N-1)
	#	push!(c, z[:,n+1+1] == tm.A*z[:,n+1] + tm.B*v[:,n+1])
	#end

	# prepare for logging if need be 
	if logging
		outfile = open("temp.csv", "w")
		save(outfile, xd)
	end

	if verbose; print_header(); end
	i = 1
	not_finished = true
	es = 0.; cs = 0.; ts = 0.; dd = 0.; step_size = 0.
	while not_finished

		# Find gradients, descent step, and step size. Then descend!
		gradients!(ad, bd, em, tm, xd, ud, start_idx)

		# make constraints
		c = Array(Constraint, 0)
		push!(c, z[:,1] == 0)
		for n = 0:(N-1)
			temp_A = eye(3)
			temp_A[1,3] = -tm.h * sin(xd[n+1][3])
			temp_A[2,3] = tm.h * cos(xd[n+1][3])
			A,B = linearize(tm, xd[n+1], ud[n+1])
			push!(c, z[:,n+1+1] == A*z[:,n+1] + B*v[:,n+1])
		end

		zd, vd = convex_descent(ad, bd, N, z, v, c, tm.R[1,1], start_idx)
		step_size = get_step_size(tm.descender, i)
		descend!(xd, ud, zd, vd, step_size, N)

		# compute statistics and report
		es, cs, ts = all_scores(em, tm, xd, ud, start_idx)
		# only temporarily commented out 2/07/2017
		#dd = directional_derivative(ad, bd, zd, vd)
		dd = 0.05
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

#function check_convergence(es::Float64, es_crit::Float64, i::Int, max_iters::Int, verbose::Bool)
#	not_finished = true
#	if es < es_crit
#		not_finished = false
#		if verbose
#			println("reached ergodic criterion...")
#		end
#	end
#	if i > max_iters
#		not_finished = false
#		if verbose
#			println("max iterations reached...")
#		end
#	end
#	return not_finished
#end
#
#function step_report(i::Int, es::Float64, cs::Float64, ts::Float64, dd::Float64, step_size::Float64)
#	@printf " %-7i" i
#	@printf " %-14.7f" es
#	@printf " %-14.7f" cs
#	@printf " %-12.7f" ts
#	@printf " %-12.7f" dd
#	@printf " %-7.5f" step_size
#	println()
#end
#
#function print_header()
#	print_dashes()
#	println(" iter  |ergodic score |control score |total score |direc deriv |step size")
#	print_dashes()
#end
#function print_dashes()
#	println("--------------------------------------------------------------------------")
#end
#
#function gradients!(ad::Matrix{Float64}, bd::Matrix{Float64}, em::ErgodicManager, tm::TrajectoryManager, xd::VVF, ud::VVF, start_idx::Int)
#	ni =  1
#	for n = 0:(tm.N-1)
#		an_x, an_y = compute_ans(em, xd, tm, n, start_idx)
#		ad[1,ni] = an_x
#		ad[2,ni] = an_y
#		# assume that only the first two state variables matter
#		# in fact, this isn't necessary if it's been initalized to zeros...
#		for i = 3:tm.n
#			ad[i,ni] = 0.0
#		end
#
#		#bd[n+1] = tm.h * tm.R * ud[n+1]
#		# TODO: allow for different tm.m
#		bd[1,ni] = tm.h * (tm.R[1,1]*ud[ni][1] + tm.R[1,2]*ud[ni][2])
#		bd[2,ni] = tm.h * (tm.R[2,1]*ud[ni][1] + tm.R[2,2]*ud[ni][2])
#
#		ni += 1
#	end
#end
#
#function compute_ans(em::ErgodicManager, xd::VVF, tm::TrajectoryManager, n::Int, start_idx::Int)
#	xnx = xd[n + start_idx + 1][1]
#	xny = xd[n + start_idx + 1][2]
#	L = em.L
#
#	an_x = 0.0
#	an_y = 0.0
#	 
#	for k1 = 0:em.K
#		for k2 = 0:em.K
#			hk = em.hk[k1+1,k2+1]
#
#			dFk_dxn1 = -k1*pi*sin(k1*pi*xnx/L)*cos(k2*pi*xny/L) / (hk*L)
#			dFk_dxn2 = -k2*pi*cos(k1*pi*xnx/L)*sin(k2*pi*xny/L) / (hk*L)
#
#			fk = 0.0
#			for i in 0:(tm.N-1)
#				x = xd[i + start_idx + 1]
#				fk += cos(k1*pi*x[1]/L) * cos(k2*pi*x[2]/L) / hk
#			end
#			c = em.Lambdak[k1+1,k2+1] * (tm.h*fk/tm.T - em.phik[k1+1,k2+1])
#			an_x += c*dFk_dxn1
#			an_y += c*dFk_dxn2
#		end
#	end
#	an_x *= 2.0*tm.h
#	an_y *= 2.0*tm.h
#	return an_x, an_y
#end
#
## descent direction using convex optimization
#function convex_descent(a::Matrix{Float64}, b::Matrix{Float64}, N::Int, z::Variable, v::Variable, c::Vector{Constraint}, r::Float64, start_idx::Int)
#
#	# create the problem
#	N_range = (0:N-1) + start_idx + 1
#	problem = minimize(vecdot(a,z[:,N_range]) + vecdot(b,v) + sumsquares(z) + r*sumsquares(v), c)
#
#	# solve the problem
#	solve!(problem, SCSSolver(verbose=0,max_iters=100000))
#
#	return z.value, v.value
#end
#
#
## modifies (xd,ud) by moving step_size in direction (zd,vd)
#function descend!(xd::VVF, ud::VVF, zd::Matrix{Float64}, vd::Matrix{Float64}, step_size::Float64, N::Int)
#	for i = 0:(N-1)
#		ud[i+1][1] += step_size*vd[1,i+1]
#		ud[i+1][2] += step_size*vd[2,i+1]
#		xd[i+1][1] += step_size*zd[1,i+1]
#		xd[i+1][2] += step_size*zd[2,i+1]
#	end
#	xd[N+1][1] += step_size*zd[1,N+1]
#	xd[N+1][2] += step_size*zd[2,N+1]
#end

# called if logging, not meant for general use
#function save(outfile::IOStream, xd::VVF)
#	n = length(xd[1])
#	for xi in xd
#		for i = 1:(n-1)
#			wi = xi[i]
#			write(outfile,"$(xi[i]),")
#		end
#		write(outfile,"$(xi[n])\n")
#	end
#end
