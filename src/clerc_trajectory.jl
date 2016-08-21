######################################################################
# clerc_trajectory.jl
# CLErC (Constrained Linear Ergodic Control)
# Assumes linear dynamics, but you can constrain norm of effort
######################################################################

using Convex, SCS

function clerc_trajectory(em::ErgodicManager, tm::TrajectoryManager; verbose=true, max_iters::Int=30)
	xd0, ud0 = initialize(tm.initializer, em, tm.x0, tm.h, tm.N)
	clerc_trajectory(em, tm, xd0, ud0; verbose=verbose,max_iters=max_iters)
end

function clerc_trajectory(em::ErgodicManager, tm::TrajectoryManager, xd0::VV_F, ud0::VV_F; verbose::Bool=true, max_iters::Int=30)

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

	# TODO: really a termination condition
	for i = 1:max_iters
		gradients!(ad, bd, em, tm, xd, ud)

		zd, vd = convex_descent(ad, bd, N, z, v, c)

		#step_size = 1.0 / i
		#step_size = 101.50 / sqrt(i)
		step_size = .15 / sqrt(i)

		# printing statistics for testing
		if verbose
			iteration_report(i,ad,bd,zd,vd,xd,ud,step_size,em,tm)
		end

		simple_project!(xd, ud, zd, vd, step_size, N)
	end
	return xd, ud
end

function iteration_report(i::Int,ad::Matrix{Float64},bd::Matrix{Float64},zd::Matrix{Float64},vd::Matrix{Float64},xd::VV_F,ud::VV_F, step_size::Float64, em::ErgodicManager, tm::TrajectoryManager)
	println("i = ",i) 
	dd = directional_derivative(ad, bd, zd, vd)
	sdd = scaled_dd(ad, bd, zd, vd)
	es = ergodic_score(em, xd)
	ts = total_score(em, xd, ud, tm.T)
	cs = control_score(ud, tm.R, tm.h)
	println("es = ", es)
	println("ts = ", ts)
	println("cs = ", cs)
	println("dd = ", dd)
	println("scaled_dd = ", sdd)
	println("alpha = ", step_size)
	println("##################################")
end

function gradients!(ad::Matrix{Float64}, bd::Matrix{Float64}, em::ErgodicManager, tm::TrajectoryManager, xd::VV_F, ud::VV_F)
	for n = (tm.N-1):-1:0
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
function convex_descent(a::Matrix{Float64}, b::Matrix{Float64}, N::Int, z::Variable, v::Variable, c::Vector{Constraint})

	# create the problem
	problem = minimize(vecdot(a,z[:,1:N]) + vecdot(b,v) + sumsquares(z) + .01*sumsquares(v), c)

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


function simple_project!(xd::VV_F, ud::VV_F, zd::Matrix{Float64}, vd::Matrix{Float64}, step_size::Float64, N::Int)
	for i = 0:(N-1)
		ud[i+1][1] += step_size*vd[1,i+1]
		ud[i+1][2] += step_size*vd[2,i+1]
		xd[i+1][1] += step_size*zd[1,i+1]
		xd[i+1][2] += step_size*zd[2,i+1]
	end
	ud[N+1] = ud[N]
end
