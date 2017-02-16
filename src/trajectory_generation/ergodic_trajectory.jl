######################################################################
# ergodic_trajectory.jl
# doesn't actually work
# currently, just a graveyard for code that might one day work
#
# clerc_trajectory has some ergodic trajectory generation, but it
#  only works for linear problems
######################################################################

######################################################################
# Creating discrete control policies
# Right now, I assume simple integrator
######################################################################
# assuming Q_n = Q
# assuming R_n = R
# assuming A_n = A
# assuming B_n = B
function gains_and_descent(em::ErgodicManager, A::Matrix{Float64}, B::Matrix{Float64}, N::Int, xd::VV_F, ud::VV_F)

	# arbitrary design parameters
	Q = 1e2*eye(2)
	R = 1e2*eye(2)
	Q = eye(2)
	R = eye(2)

	# avoid repeated 
	AT = A'
	BT = B'

	# required for gains and descent direction
	K = Array(Matrix{Float64}, N)
	P = Array(Matrix{Float64}, N+1)
	P[N+1] = Q

	# required for descent direction
	r = Array(Vector{Float64}, N+1)
	r[N+1] = compute_an(em, xd, N, N)
	v = Array(Vector{Float64}, N+1)
	z = Array(Vector{Float64}, N+1)

	# partial derivatives of objective function J
	# TODO: what is the right length for these guys???
	a = Array(Vector{Float64}, N)
	b = Array(Vector{Float64}, N)

	for n = (N-1):-1:0

		# required for gains and descent direction
		G = R + BT*P[n+1+1]*B
		Kn = inv(G)*BT*P[n+1+1]*A
		#Kn = -inv(G)*BT*P[n+1+1]*A		# different from paper but right
		KnT = Kn'
		K[n+1] = Kn
		P[n+1] = Q + AT*P[n+1+1]*A - KnT*G*Kn

		# required for descent direction
		a[n+1] = compute_an(em, xd, N, n)
		b[n+1] = compute_bn(ud[n+1])
		r[n+1] = (AT - KnT*BT)*r[n+1+1] + a[n+1] - KnT*b[n+1]
	end

	# compute descent direction for position and velocity
	# now that we've computed r, we need to compute v and z
	z[1] = zeros(2)
	for n = 0:(N-1)
		v[n+1] = b[n+1] + BT*P[n+1]*z[n+1] + BT*r[n+1+1]
		z[n+1+1] = A*z[n+1] + B*v[n+1]
	end
	v[N+1] = v[N]		# last control input doesn't matter (no effect)

	#return K, z, v, a, b
	return K, z, v, a, b, P
end

# like the above, but uses convex for descent direction
function gains_and_descent2(em::ErgodicManager, A::Matrix{Float64}, B::Matrix{Float64}, N::Int, xd::VV_F, ud::VV_F)

	# arbitrary design parameters
	Q = eye(2)
	R = eye(2)

	# avoid repeated 
	AT = A'
	BT = B'

	# required for gains and descent direction
	K = Array(Matrix{Float64}, N)
	P = Array(Matrix{Float64}, N+1)
	P[N+1] = Q

	# required for descent direction
	r = Array(Vector{Float64}, N+1)
	r[N+1] = compute_an(em, xd, N, N)
	v = Array(Vector{Float64}, N+1)
	z = Array(Vector{Float64}, N+1)

	# partial derivatives of objective function J
	# TODO: what is the right length for these guys???
	a = Array(Vector{Float64}, N)
	b = Array(Vector{Float64}, N)

	for n = (N-1):-1:0

		# required for gains and descent direction
		G = R + BT*P[n+1+1]*B
		Kn = inv(G)*BT*P[n+1+1]*A
		#Kn = -inv(G)*BT*P[n+1+1]*A		# different from paper but right
		KnT = Kn'
		K[n+1] = Kn
		P[n+1] = Q + AT*P[n+1+1]*A - KnT*G*Kn

		# required for descent direction
		a[n+1] = compute_an(em, xd, N, n)
		b[n+1] = compute_bn(ud[n+1])
		r[n+1] = (AT - KnT*BT)*r[n+1+1] + a[n+1] - KnT*b[n+1]
	end

	z, v = convex_descent(A, B, a, b, N)

	#return K, z, v, a, b
	return K, z, v, a, b, P
end

# descent direction using convex optimization
function convex_descent(A::Matrix{Float64}, B::Matrix{Float64}, av::VV_F, bv::VV_F, N::Int)
	# convert av and bv to better format (array instead of vec of vecs)
	a = zeros(2,N)
	b = zeros(2,N)
	for i = 1:N
		a[:,i] = av[i]
		b[:,i] = bv[i]
	end

	# variables
	z = Variable(2, N+1)
	v = Variable(2, N)

	# constraints
	c = Array(Constraint, 0)
	push!(c, z[:,1] == 0)
	for n = 0:(N-1)
		push!(c, z[:,n+1+1] == A*z[:,n+1] + B*v[:,n+1])
	end

	# create the problem
	problem = minimize(vecdot(a,z[:,1:N]) + vecdot(b,v) + sumsquares(z) + .01*sumsquares(v), c)

	# solve the problem
	solve!(problem, SCSSolver(verbose=0,max_iters=100000))

	# put z and v into vec of vec format
	vd = Array(Vector{Float64}, N+1)
	zd = Array(Vector{Float64}, N+1)
	for i = 1:N
		zd[i] = [z.value[1,i], z.value[2,i]]
		vd[i] = [v.value[1,i], v.value[2,i]]
	end
	# just repeat last control, it technically doesn't matter
	zd[N+1] = [z.value[1,N+1], z.value[2,N+1]]
	vd[N+1] = [v.value[1,N], v.value[2,N]]

	return zd, vd
end

function compute_bn(ui::Vector{Float64})
	h = 0.5
	R = 0.01 * eye(2)
	return h * R * ui
end

function make_trajectory(em::ErgodicManager, tm::TrajectoryManager; verbose=true, max_iters::Int=30)
	xd0, ud0 = initialize(tm.initializer, em, tm.x0, tm.h, tm.N)
	make_trajectory(em, xd0, ud0, tm.h; verbose=verbose,max_iters=max_iters)
end

function make_trajectory(em::ErgodicManager, xd::VV_F, ud::VV_F, h::Float64; verbose::Bool=true, max_iters::Int=30)

	N = length(xd) - 1
	T = N*h

	# generate linearized dynamics (assumed constant)
	# TODO: generate these from system dynamics
	A = eye(2)
	B = h*eye(2)

	# TODO: really a termination condition
	for i = 1:max_iters
		K, zd, vd, ad, bd, P = gains_and_descent2(em, A, B, N, xd, ud)
		#if i > 2
		#	step_size = armijo_ls(em, xd, ud, zd, vd, ad, bd, T)
		#	#step_size = 0.01
		#else
		#	step_size = 1.0
		#	#step_size = armijo_ls(em, xd, ud, zd, vd, ad, bd, T)
		#end
		#step_size = armijo_ls(em, xd, ud, zd, vd, ad, bd, T)
		#step_size = 1.0
		#step_size = 1.0 / i
		#step_size = 101.50 / sqrt(i)
		step_size = .15 / sqrt(i)

		# printing statistics for testing
		if verbose
			println("i = ",i) 
			dd = directional_derivative(ad, bd, zd, vd)
			sdd = scaled_dd(ad, bd, zd, vd)
			es = ergodic_score(em, xd)
			ts = total_score(em, xd, ud, T)
			cs = control_score(ud, 0.01*eye(2), h)
			println("es = ", es)
			println("ts = ", ts)
			println("cs = ", cs)
			println("dd = ", dd)
			println("scaled_dd = ", sdd)
			println("alpha = ", step_size)
			println("##################################")
		end

		#project!(xd, ud, zd, vd, step_size, K, N, h)
		project2!(xd, ud, zd, vd, step_size, K, N, h)
	end
	return xd, ud
	#return xd, ud, zd, vd, K, ad, bd, P
end

# TODO: determine if m should involve scaled directional derivative
function armijo_ls(em, xd, ud, zd, vd, ad, bd, T)
	print("armijo starting...")
	tau = 0.5
	c = 0.5
	#c = 0.9
	alpha = 1e4

	# compute m = p' * grad f(x)
	#m = directional_derivative(ad, bd, zd, vd)
	normalizer!(zd, vd)
	m = directional_derivative(ad, bd, zd, vd)
	#m = scaled_dd(ad, bd, zd, vd)

	f_x = total_score(em, xd, ud, T)

	while total_score(em, xd, ud, zd, vd, alpha, T) > f_x + alpha*c*m
		alpha *= tau
	end
	println("done.")
	return alpha
end


# Modifies trajectory with step size and descent direction.
# Then projects this modified trajectory to a feasible trajectory.
# TODO: is this necessary for linear dynamics?
#
# Process:
#  alpha_n = xd_n + step_size*zd_n
#  mu_n    = ud_n + step_size*vd_n
#
#  ud_n = mu_n + K(alpha_n - xd_n)
#
#  substituting leads to:
#
#  ud_n += step_size*vd_n + K(step_size*zd_n)
#
function project!(xd::VV_F, ud::VV_F, zd::VV_F, vd::VV_F, step_size::Float64, K::VMF64, N::Int, h::Float64)
	for i = 0:(N-1)
		ud[i+1] += step_size*vd[i+1] + K[i+1]*step_size*zd[i+1]
		xd[i+1+1] = xd[i+1] + h*ud[i+1]

		# Keep it in bounds (do we really have to do this?
		if false
			if xd[i+1+1][1] < 0.0
				xd[i+1+1][1] = 0.0
				ud[i+1][1] = (xd[i+1+1][1] - xd[i+1][1])/h
			end
			if xd[i+1+1][2] < 0.0
				xd[i+1+1][2] = 0.0
				ud[i+1][2] = (xd[i+1+1][2] - xd[i+1][2])/h
			end
			if xd[i+1+1][1] > 1.0
				xd[i+1+1][1] = 1.0
				ud[i+1][1] = (xd[i+1+1][1] - xd[i+1][1])/h
			end
			if xd[i+1+1][2] > 1.0
				xd[i+1+1][2] = 1.0
				ud[i+1][2] = (xd[i+1+1][2] - xd[i+1][2])/h
			end
		end
	end
	ud[N+1] = ud[N]
end

# Doesn't project, just modifies
# for testing
function project2!(xd::VV_F, ud::VV_F, zd::VV_F, vd::VV_F, step_size::Float64, K::VMF64, N::Int, h::Float64)
	for i = 0:(N-1)
		ud[i+1] += step_size*vd[i+1]
		xd[i+1] += step_size*zd[i+1]

	end
	ud[N+1] = ud[N]
end

# X is a trajectory
function compute_an(em::ErgodicManager, xd::VV_F, N::Int, n::Int)
	xnx = xd[n+1][1]
	xny = xd[n+1][2]
	L = em.L
	an = zeros(2)

	h = 0.5	  # TODO: input parameter
	T = h*N	  # TODO: input parameter (is this just h*N?)
	 
	for k1 = 0:em.K
		for k2 = 0:em.K
			hk = em.hk[k1+1,k2+1]

			dFk_dxn1 = -k1*pi*sin(k1*pi*xnx/L)*cos(k2*pi*xny/L) / (hk*L)
			dFk_dxn2 = -k2*pi*cos(k1*pi*xnx/L)*sin(k2*pi*xny/L) / (hk*L)

			fk = 0.0
			for i = 0:(N-1)
				x = xd[i+1]
				fk += cos(k1*pi*x[1]/L) * cos(k2*pi*x[2]/L) / hk
			end
			c = em.Lambdak[k1+1,k2+1] * (h*fk/T - em.phik[k1+1,k2+1])
			an[1] += c*dFk_dxn1
			an[2] += c*dFk_dxn2
		end
	end
	an[1] *= 2.0*h
	an[2] *= 2.0*h
	return an
end
