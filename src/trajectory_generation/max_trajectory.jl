######################################################################
# max_trajectory.jl
# tries to maximize phi(x) collected
# TODO: this is a wreck. clean it so it looks like clerc_trajectory.
######################################################################

function max_trajectory(tm::TrajectoryManager, mu::VF, Sigma::MF)

	# initialize trajectory
	# TODO do this correctly
	xd, ud = initialize_trajectory(N, h, x0)
	xd0, ud0 = initialize(tm.initializer, em, tm)

	# generate linearized dynamics (assumed constant)

	# TODO: really a termination condition
	for i = 1:100
		ad, bd = compute_gradients(tm, xd, ud, mu, Sigma)
		zd, vd = convex_descent(A, B, ad, bd, N)
		step_size = .15 / sqrt(i)

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

		project3!(xd, ud, zd, vd, step_size)
	end

	return xd, ud
end

export max_trajectory


# TODO: should really take in a bunch of gaussians here
function compute_gradients(tm::TrajectoryManager, xd::VVF, ud::VVF, mu::VF, Sigma::MF)
	N = length(xd) - 1

	a = Array(VF, N+1)
	b = Array(VF, N)

	for ni = 1:N
		a[ni] = compute_an(mu, Sigma, xd[ni], tm.h)
		b[ni] = tm.h * tm.R * ud[ni]
	end
	a[N+1] = compute_an(mu, Sigma, xd[N+1], tm.h)

	return a, b
end

function compute_an(mu::VF, Sigma::MF, x::VF, h::Float64)
	xmu = x - mu
	inv_Sigma = inv(Sigma)
	q = 1.0
	c = q*h*exp(-.5dot(xmu, inv_Sigma*xmu)) / (2*pi*sqrt(det(Sigma)))
	return c*inv_Sigma*xmu
end



function project3!(xd::VVF, ud::VVF, zd::VVF,vd::VVF,step_size::Float64)
	N = length(ud) - 1
	for i = 0:(N-1)
		ud[i+1] += step_size*vd[i+1]
		xd[i+1] += step_size*zd[i+1]
	end
	ud[N+1] = ud[N]
end
