######################################################################
# new_trajectory.jl
#
# Here I try to apply the original way of doing it.
######################################################################

export new_trajectory

function new_trajectory(em::ErgodicManager, tm::TrajectoryManager; verbose::Bool=true, logging::Bool=false, max_iters::Int=30, es_crit::Float64=0.003,right::Bool=true)
	xd0, ud0 = initialize(tm.initializer, em, tm)
	new_trajectory(em, tm, xd0, ud0; verbose=verbose, logging=logging, max_iters=max_iters, es_crit=es_crit,right=right)
end

function new_trajectory(em::ErgodicManager, tm::TrajectoryManager, xd0::VV_F, ud0::VV_F; verbose::Bool=true, logging::Bool=false, max_iters::Int=30, es_crit::Float64=0.003, right::Bool=true)

	# let's not overwrite the initial trajectories
	xd = deepcopy(xd0)
	ud = deepcopy(ud0)
	N = tm.N
	start_idx = right ? 1 : 0

	# matrices for gradients
	ad = zeros(tm.dynamics.n, N)
	bd = zeros(tm.dynamics.m, N)

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

		# determine gradients used in optimization
		gradients!(ad, bd, em, tm, xd, ud, start_idx)

		# Find gains K and descent direction (LQ)
		A,B = linearize(tm.dynamics, xd, ud, tm.h)
		#A = tm.dynamics.A
		#B = tm.dynamics.B
		#K, C = LQ(tm.dynamics.A, tm.dynamics.B, ad, bd, tm.Q, tm.R, tm.N)
		K, C = LQ(A, B, ad, bd, tm.Q, tm.R, tm.N)
		#zd, vd = apply_LQ_gains(tm.dynamics.A, tm.dynamics.B, K, C)
		zd, vd = apply_LQ_gains(A, B, K, C)

		# determine step size and descend
		step_size = get_step_size(tm.descender, em, tm, xd, ud, zd, vd, ad, bd, K, i)

		# descend and project
		xd, ud = project(em, tm, K, xd, ud, zd, vd, step_size)
		#descend!(xd, ud, zd, vd, step_size, N)

		# compute statistics and report
		es, cs, ts = all_scores(em, tm, xd, ud, start_idx)
		dd = directional_derivative(ad, bd, zd, vd)
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
