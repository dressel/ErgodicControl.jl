######################################################################
# max_trajectory.jl
# tries to maximize phi(x) collected
# TODO: this is a wreck. clean it so it looks like clerc_trajectory.
######################################################################

function max_trajectory(em::ErgodicManager, tm::TrajectoryManager, mu::VF, Sigma::MF)
	return max_trajectory(em, tm, [mu], [Sigma])
end

function max_trajectory(em::ErgodicManager, tm::TrajectoryManager, mus::VVF, Sigmas::VMF, weights::VF=ones(length(mus)); max_iters=100)

	# initialize trajectory
	xd, ud = initialize(tm.initializer, em, tm)

	# TODO: really a termination condition
	for i = 1:max_iters

		ad, bd = compute_gradients(tm, xd, ud, mus, Sigmas, weights)
		A, B = linearize(tm.dynamics, xd, ud, tm.h)
		K, C = LQ(A, B, ad, bd, tm.Qn, tm.Rn, tm.N)

		#zd, vd = convex_descent(A, B, ad, bd, N)
		zd, vd = apply_LQ_gains(A, B, K, C)
		step_size = .15 / sqrt(i)
		#step_size = .25 / sqrt(i)

		# printing statistics for testing
		println("i = ",i) 
		dd = directional_derivative(ad, bd, zd, vd)
		sdd = scaled_dd(ad, bd, zd, vd)
		#es = ergodic_score(em, xd)
		#ts = total_score(em, xd, ud, T)
		#cs = control_score(ud, 0.01*eye(2), T)
		#println("es = ", es)
		#println("ts = ", ts)
		#println("cs = ", cs)
		println("dd = ", dd)
		println("scaled_dd = ", sdd)
		#println("alpha = ", step_size)
		println("##################################")

		descend!(xd, ud, zd, vd, step_size)
	end

	return xd, ud
end

export max_trajectory


# TODO: should really take in a bunch of gaussians here
# returns a and b, each of which is an array of vectors
function compute_gradients(tm::TrajectoryManager, xd::VVF, ud::VVF, mus::VVF, Sigmas::VMF, weights::VF)
	#a = Array(VF, tm.N+1)
	#b = Array(VF, tm.N)
	a = zeros(tm.dynamics.n, tm.N+1)
	b = zeros(tm.dynamics.m, tm.N)

	for ni = 1:tm.N
		a[1:2, ni] = compute_an(mus, Sigmas, weights, xd[ni], tm.h, tm.q)
		b[:, ni] = tm.h * tm.R * ud[ni]
	end
	a[1:2,tm.N+1] = compute_an(mus, Sigmas, weights, xd[tm.N+1], tm.h, tm.q)

	return a, b
end

function compute_an(mus::VVF, Sigmas::VMF, weights::VF, x::VF, h::Float64, q::Float64)
	num_gauss = length(weights)
	an = zeros(2)
	for i = 1:num_gauss
		mu = mus[i]
		Sigma = Sigmas[i]
		xmu = x[1:2] - mu
		inv_Sigma = inv(Sigma)
		c = q*h*exp(-.5dot(xmu, inv_Sigma*xmu)) / (2*pi*sqrt(det(Sigma)))
		an += weights[i]*c*inv_Sigma*xmu
	end
	return an
end
