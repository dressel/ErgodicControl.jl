######################################################################
# dynamics.jl
# Each trajectory manager will have a type of dynamics
######################################################################
export Dynamics, LinearDynamics, DubinsDynamics, linearize

type LinearDynamics <: Dynamics
	n::Int
	m::Int
	A::Matrix{Float64}
	B::Matrix{Float64}

	function LinearDynamics(A, B)
		n,m = size(B)
		return new(n,m,deepcopy(A),deepcopy(B))
	end
end

type DubinsDynamics <: Dynamics
	n::Int
	m::Int

	v0::Float64		# constant speed
	r::Float64		# minimum turn radius

	# Constructors
	DubinsDynamics() = DubinsDynamics(1.0, 1.0)
	DubinsDynamics(v0::Real, r::Real) = new(3, 1, v0, r)
end


function linearize(d::Dynamics, x::VVF, u::VVF, h::Float64)
	N = length(u)
	A = VMF()
	B = VMF()
	for n = 1:N
		An, Bn = linearize(d, x[n], u[n], h)
		push!(A, An)
		push!(B, Bn)
	end
	return A,B
end


function linearize(ld::LinearDynamics, x::VF, u::VF, h::Float64)
	return ld.A, ld.B
end

function linearize(ld::DubinsDynamics, x::VF, u::VF, h::Float64)
	A = eye(3)
	A[1,3] = -h * sin(x[3]) * ld.v0
	A[2,3] = h * cos(x[3]) * ld.v0

	B = zeros(3,1)
	B[3] = h/ld.r

	return A, B
end


######################################################################
# forward_euler
######################################################################
function forward_euler(tm::TrajectoryManager, x::VF, u::VF)
	forward_euler(tm.dynamics, x, u, tm.h)
end

function forward_euler(tm::TrajectoryManager, ud::VVF)
	N = length(ud)
	xd = Array(Vector{Float64}, N+1)

	xd[1] = deepcopy(tm.x0)
	for i = 1:tm.N
		xd[i+1] = forward_euler(tm.dynamics, x, ud[i], tm.h)
	end

	return xd
end

function forward_euler(ld::LinearDynamics, x::VF, u::VF, h::Float64)
	return ld.A*x + ld.B*u
end
function forward_euler(dd::DubinsDynamics, x::VF, u::VF, h::Float64)
	xp = deepcopy(x)
	xp[1] += cos(x[3]) * dd.v0 * h
	xp[2] += sin(x[3]) * dd.v0 * h
	u_val = u[1]
	if u_val > 1.
		u_val = 1.
	end
	if u_val < -1.
		u_val = -1.
	end
	xp[3] += u_val / dd.r * h
	return xp
end
