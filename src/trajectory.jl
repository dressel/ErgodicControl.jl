######################################################################
# trajectory.jl
#
# handles the generation of ergodic trajectories
######################################################################

using Convex, SCS

type TrajectoryManager
	N::Int
	h::Float64
	T::Float64
	Q::Matrix{Float64}
	R::Matrix{Float64}
	q::Float64

	function TrajectoryManager(N::Int, h::Float64)
		tm = new()
		tm.N = N
		tm.h = h
		tm.T = N*h
		tm.Q = eye(2)
		tm.R = 0.01 * eye(2)
		tm.q = 1.0
		return tm
	end
end


"""
`decompose(em, traj::VV_F)`

`decompose(em, traj::V_T2F)`

Decomposes a set of positions into a set of `ck` Fourier coefficients.
"""
function decompose(em::ErgodicManager, traj::VV_F)
	traj2 = [(traj[i][1], traj[i][2]) for i = 1:length(traj)]
	return decompose(em, traj2)
end
function decompose(em::ErgodicManager, traj::V_T2F)
	K = em.K
	N = length(traj)-1
	ck = zeros(K+1, K+1)
	for k1 = 0:K
		kpiL1 = k1 * pi / em.L
		for k2 = 0:K
			kpiL2 = k2 * pi / em.L
			hk = em.hk[k1+1, k2+1]
			fk_sum = 0.0
			# now loop over time
			for n = 0:N-1
				xn = traj[n+1]
				fk_sum += cos(kpiL1 * xn[1])  * cos(kpiL2 * xn[2])
			end
			ck[k1+1, k2+1] = fk_sum / (hk * N)
		end
	end
	return ck
end


"""
`ergodic_score(em, traj::V_T2F)`

First breaks down the trajectory into components ck.
"""
function ergodic_score(em::ErgodicManager, traj::V_T2F)
	ck = decompose(em, traj)
	return ergodic_score(em, ck)
end
function ergodic_score(em::ErgodicManager, traj::VV_F)
	ck = decompose(em, traj)
	return ergodic_score(em, ck)
end
function ergodic_score(em::ErgodicManager, ck::Matrix{Float64})
	val = 0.0
	for k1 = 0:em.K
		for k2 = 0:em.K
			d = em.phik[k1+1,k2+1] - ck[k1+1,k2+1]
			val += em.Lambdak[k1+1,k2+1] * d * d
		end
	end
	return val
end

"""
`control_score(ud::VV_F, R, T)`

Assumes only non-zero elements of `R` are corners.
"""
function control_score(ud::VV_F, R::Matrix{Float64}, T::Float64)
	N = length(ud) - 1
	h = T / N
	cs = 0.0
	for ui in ud[1:end-1]
		cs += R[1,1] * ui[1] * ui[1]
		cs += R[2,2] * ui[2] * ui[2]
	end
	return 0.5 * h * cs
end


"""
`total_score(em, xd::VV_F, ud::VV_F, T::Float64)`

Computes the total score (q*ergodic_score + sum_n h/2 un'Rn un)
"""
# TODO: actually get q and R from the correct place 
function total_score(em::ErgodicManager, xd::VV_F, ud::VV_F, T::Float64)
	q = 1.0
	R = 0.01 * eye(2)
	return q * ergodic_score(em, xd) + control_score(ud, R, T)
end
# TODO: let's not make this so shitty...
function total_score(em::ErgodicManager, xd::VV_F, ud::VV_F, zd::VV_F, vd::VV_F, alpha::Float64, T::Float64)
	xd2 = deepcopy(xd)
	ud2 = deepcopy(ud)
	for i = 1:length(xd2)
		xd2[i][1] += alpha * zd[i][1]
		xd2[i][2] += alpha * zd[i][2]

		ud2[i][1] += alpha * vd[i][1]
		ud2[i][2] += alpha * vd[i][2]
	end
	return total_score(em, xd2, ud2, T)
end


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
			#for i = 1:(N-1)
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

function compute_bn(ui::Vector{Float64})
	h = 0.5
	R = 0.01 * eye(2)
	return h * R * ui
end

# unscaled
# g_f1 = grad_1 f
# g_f2 = grad_2 f
# u1, u2 are components of direction
function directional_derivative(g_f1::VV_F, g_f2::VV_F, u1::VV_F, u2::VV_F)
	N = length(g_f1)
	dd = 0.0
	for i = 1:N
		dd += dot(g_f1[i], u1[i]) + dot(g_f2[i], u2[i])
	end
	return dd
end
# like above, but scales the direction first
function scaled_dd(g_f1::VV_F, g_f2::VV_F, u1::VV_F, u2::VV_F)
	dd = directional_derivative(g_f1, g_f2, u1, u2)
	norm_factor = sqrt( dot(u1,u1) + dot(u2,u2) )
	return dd / norm_factor
end

function make_trajectory(em::ErgodicManager, N::Int, h::Float64, x0::T2F)
	xd, ud = initialize_trajectory(N, h, x0)
	make_trajectory(em, xd, ud, h)
end
function make_trajectory(em::ErgodicManager, xd::VV_F, ud::VV_F, h::Float64)

	N = length(xd) - 1
	T = N*h

	# initialize trajectory
	#xd, ud = initialize_trajectory(N, h, x0)

	# generate linearized dynamics (assumed constant)
	# TODO: generate these from system dynamics
	A = eye(2)
	B = h*eye(2)

	# TODO: really a termination condition
	zd = 3
	vd = 3
	K = 3
	ad = 3
	bd = 3
	P = 3
	for i = 1:1000
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
		println("i = ",i) 
		dd = directional_derivative(ad, bd, zd, vd)
		sdd = scaled_dd(ad, bd, zd, vd)
		es = ergodic_score(em, xd)
		ts = total_score(em, xd, ud, T)
		cs = control_score(ud, 0.01*eye(2), T)
		println("es = ", es)
		println("ts = ", ts)
		println("cs = ", cs)
		println("dd = ", dd)
		println("scaled_dd = ", sdd)
		println("alpha = ", step_size)
		println("##################################")

		#project!(xd, ud, zd, vd, step_size, K, N, h)
		project2!(xd, ud, zd, vd, step_size, K, N, h)
	end
	#return xd, ud
	return xd, ud, zd, vd, K, ad, bd, P
end

# normalizes zd and vd
# recall that zd and vd are together a direction, we normalize em both
function normalizer!(zd::VV_F, vd::VV_F)
	norm_factor = sqrt(dot(zd,zd) + dot(vd,vd))
	for i = 1:length(zd)
		zd[i][1] /= norm_factor
		zd[i][2] /= norm_factor
		vd[i][1] /= norm_factor
		vd[i][2] /= norm_factor
	end
end


# currently, I just move in a random direction
# perhaps find direction to mean and move towards it
function initialize_trajectory(N::Int, h::Float64, x0::T2F)
	xd = [[x0[1], x0[2]] for i = 1:N+1]
	#ud = [.01*ones(2) for i = 1:N+1]
	#ud = [[.01,0.] for i = 1:N+1]
	ud = [[.01,0.01] for i = 1:N+1]
	for i = 2:N+1
		xd[i][1] = xd[i-1][1] + h*ud[i-1][1]
		xd[i][2] = xd[i-1][2] + h*ud[i-1][2]
	end
	return xd, ud
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
		#xd[i+1+1] = xd[i+1] + h*ud[i+1]
		xd[i+1] += step_size*zd[i+1]

		# Keep it in bounds (do we really have to do this?
		#if xd[i+1+1][1] < 0.0
		#	xd[i+1+1][1] = 0.0
		#	ud[i+1][1] = (xd[i+1+1][1] - xd[i+1][1])/h
		#end
		#if xd[i+1+1][2] < 0.0
		#	xd[i+1+1][2] = 0.0
		#	ud[i+1][2] = (xd[i+1+1][2] - xd[i+1][2])/h
		#end
		#if xd[i+1+1][1] > 1.0
		#	xd[i+1+1][1] = 1.0
		#	ud[i+1][1] = (xd[i+1+1][1] - xd[i+1][1])/h
		#end
		#if xd[i+1+1][2] > 1.0
		#	xd[i+1+1][2] = 1.0
		#	ud[i+1][2] = (xd[i+1+1][2] - xd[i+1][2])/h
		#end
	end
	ud[N+1] = ud[N]
end

"""
`collect_info(em::ErgodicManager, traj::VV_F, d_rate::Float64)`

`collect_info(em::ErgodicManager, traj::VV_F)`

modifies em.phi according to some submodular.
we don't use the last point in the trajectory

decreases at rate `D/T`
If you spend `h` time there, it is equivalent to `h*D/(h*N) = D/N`

Returns total info picked up (a scalar value).
"""
function collect_info(em::ErgodicManager, traj::VV_F)
	N = length(traj) - 1
	D = sum(em.phi)
	d_rate = D/N
	collect_info(em, traj, d_rate)
end

function collect_info(em::ErgodicManager, traj::VV_F, d_rate::Float64)
	N = length(traj) - 1
	total_info = 0.0
	for n = 0:(N-1)
		xi,yi = find_cell(em, traj[n+1])

		# if there is enough info, grab that shit yo
		info_value = min(em.phi[xi,yi], d_rate)
		em.phi[xi,yi] -= info_value
		total_info += info_value
	end
	return total_info
end
export collect_info

function find_cell(em::ErgodicManager, x::Vector{Float64})
	x1 = round(Int, x[1] / em.cell_size, RoundDown) + 1
	x2 = round(Int, x[2] / em.cell_size, RoundDown) + 1
	if x1 > em.bins; x1 -= 1; end
	if x2 > em.bins; x2 -= 1; end
	return x1, x2
end
