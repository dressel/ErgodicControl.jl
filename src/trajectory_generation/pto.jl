######################################################################
# pto_trajectory.jl
#
# Here I try to apply the original way of doing it.
######################################################################

export pto_trajectory

function pto_trajectory(em::ErgodicManager, tm::TrajectoryManager; verbose::Bool=true, logging::Bool=false, max_iters::Int=100, es_crit::Float64=0.0, dd_crit::Float64=1e-6)
	xd0, ud0 = initialize(tm.initializer, em, tm)
	pto_trajectory(em, tm, xd0, ud0; verbose=verbose, logging=logging, max_iters=max_iters, es_crit=es_crit, dd_crit = dd_crit)
end

function pto_trajectory(em::ErgodicManager, tm::TrajectoryManager, xd0::VVF, ud0::VVF; verbose::Bool=true, logging::Bool=false, max_iters::Int=100, es_crit::Float64=0.0, dd_crit::Float64=1e-6)

	# let's not overwrite the initial trajectories
	xd = deepcopy(xd0)
	ud = deepcopy(ud0)
	N = tm.N

	# matrices for gradients
	ad = zeros(tm.dynamics.n, N+1)
	bd = zeros(tm.dynamics.m, N)

	# prepare for logging if need be 
	if logging
		outfile = open("temp.csv", "w")
		save(outfile, xd)
	end

	if verbose; print_header(); end
	i = 1
	not_finished = true
	es = 0.; cs = 0.; ts = 0.; dd = 0.; step_size = 0.; es_diff=0.0
	es_prev = 0.0
	es_count = 1
	while not_finished

		# determine gradients used in optimization
		gradients!(ad, bd, em, tm, xd, ud)

		# find gains K and descent direction (zd, vd) using LQ
		A, B = linearize(tm.dynamics, xd, ud, tm.h)
		K, C = LQ(A, B, ad, bd, tm.Qn, tm.Rn, tm.N)
		zd, vd = apply_LQ_gains(A, B, K, C)

		# determine step size and descend
		step_size = get_step_size(tm.descender, em, tm, xd, ud, zd, vd, ad, bd, K, i)

		# descend and project
		xd, ud = project(em, tm, K, xd, ud, zd, vd, step_size)

		# compute statistics and report
		es, cs, ts = all_scores(em, tm, xd, ud)
		dd = directional_derivative(ad, bd, zd, vd)
		if verbose; step_report(i, es, cs, ts, dd, step_size); end
		if logging; save(outfile, xd); end

		# check convergence
		i += 1
		es_count = abs(es - es_prev) < 1e-7 ? es_count + 1 : 0
		es_prev = es
		not_finished = check_convergence(es,es_crit,i,max_iters,dd,dd_crit,verbose, es_count)
	end

	# now that we are done, print a special finished report
	if verbose
		print_header()
		if verbose; step_report(i-1, es, cs, ts, dd, step_size); end
	end

	if logging; close(outfile); end

	return xd, ud
end
