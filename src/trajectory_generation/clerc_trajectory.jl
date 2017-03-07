######################################################################
# clerc_trajectory.jl
# CLErC (Constrained Linear Ergodic Control)
# Assumes linear dynamics, but you can constrain norm of effort
######################################################################

using Convex, SCS

function clerc_trajectory(em::ErgodicManager, tm::TrajectoryManager; verbose::Bool=true, logging::Bool=false, max_iters::Int=30, es_crit::Float64=0.003,right::Bool=true)
	xd0, ud0 = initialize(tm.initializer, em, tm)
	clerc_trajectory(em, tm, xd0, ud0; verbose=verbose, logging=logging, max_iters=max_iters, es_crit=es_crit,right=right)
end

function clerc_trajectory(em::ErgodicManager, tm::TrajectoryManager, xd0::VVF, ud0::VVF; verbose::Bool=true, logging::Bool=false, max_iters::Int=30, es_crit::Float64=0.003, right::Bool=true)

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
	for n = 0:(N-1)
		push!(c, z[:,n+1+1] == tm.dynamics.A*z[:,n+1] + tm.dynamics.B*v[:,n+1])
	end

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
		zd, vd = convex_descent(ad, bd, N, z, v, c, tm.R[1,1], start_idx)
		step_size = .15 / sqrt(i)
		#step_size = .01 / sqrt(i)
		descend!(xd, ud, zd, vd, step_size, N)

		# compute statistics and report
		es, cs, ts = all_scores(em, tm, xd, ud, start_idx)
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




# descent direction using convex optimization
function convex_descent(a::Matrix{Float64}, b::Matrix{Float64}, N::Int, z::Variable, v::Variable, c::Vector{Constraint}, r::Float64, start_idx::Int)

	# create the problem
	N_range = (0:N-1) + start_idx + 1
	problem = minimize(vecdot(a,z[:,N_range]) + vecdot(b,v) + sumsquares(z) + r*sumsquares(v), c)

	# solve the problem
	solve!(problem, SCSSolver(verbose=0,max_iters=100000))

	return z.value, v.value
end
