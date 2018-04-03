######################################################################
# dynamics.jl
# Each trajectory manager will have a type of dynamics
######################################################################
export Dynamics, LinearDynamics, DubinsDynamics, linearize
export IntegrationScheme, ForwardEuler, SymplecticEuler
export forward_euler
export GroupDynamics
export vvf2vvvf

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

type GroupDynamics <: Dynamics
	n::Int
	m::Int
	num_agents::Int
	array

	function GroupDynamics(array)
		num_agents = length(array)
		n = m = 0
		for j = 1:num_agents
			n += array[j].n
			m += array[j].m
		end
		return new(n, m, num_agents, array)
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

# Not really types, 
export SingleIntegrator, DoubleIntegrator
function SingleIntegrator(n::Int, h::Float64)
	A = eye(n)
	B = h*eye(n)
	return LinearDynamics(A, B)
end
function DoubleIntegrator(n::Int, h::Float64)
	A = eye(2*n)
	for i = 1:n
		A[i, n+i] = h
	end

	B = zeros(2*n, n)
	for i = 1:n
		B[n+i, i] = h
	end

	return LinearDynamics(A, B)
end


"""
`dynamics!(tm::TrajectoryManager, d::Dynamics)`

Not only sets `tm.dynamics = d`, but also sets reward matrices `tm.Qn`, `tm.R`, and `tm.Rn` to default matrices of the correct sizes. These defaults are:

`tm.Qn = eye(d.n)`

`tm.R = 0.01 * eye(d.m)`

`tm.Rn = eye(d.m)`

"""
function dynamics!(tm::TrajectoryManager, d::Dynamics)
	tm.dynamics = d
	tm.Qn = eye(d.n)
	tm.R = 0.01 * eye(d.m)
	tm.Rn = eye(d.m)
end


######################################################################
# linearization
######################################################################

# linearizes about a trajectory
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

#function linearize(vtm::Vector{TrajectoryManager}, xds::VVVF, uds::VVVF)
#	num_agents = length(vtm)
#	N = length(xds[1]) - 1
#
#	# Make A and B a matrix with N rows and num_agents cols
#	A, B = linearize(vmt[1].dynamics, xds[1], uds[1], vmt[1].h)
#	for j = 2:num_agents
#		tm = vtm[j]
#		At, Bt = linearize(vtm[j].dynamics, xds[j], uds[j], vtm[j].h)
#		A = hcat(A, At)
#		B = hcat(B, Bt)
#	end
#
#	# now concatenate along each of the rows
#	dims = ones(Int, num_agents)
#	dims[1] = 1
#	for n = 1:N
#		push!(bigA, cat(dims, A[n]...))
#		push!(bigB, cat(dims, B[n]...))
#	end
#
#	return bigA, bigB
#end


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

function gsplit(gd::GroupDynamics, x::VF)
	xarr = VVF()
	xind = 1
	for j = 1:gd.num_agents
		d = gd.array[j]
		push!(xarr, x[xind + d.n - 1])
		xind += d.n
	end
	return xarr
end

function gsplit(gd::GroupDynamics, x::VF, u::VF)
	xarr = VVF()
	uarr = VVF()
	xind = uind = 1
	for j = 1:gd.num_agents
		d = gd.array[j]
		push!(xarr, x[xind:(xind + d.n - 1)])
		push!(uarr, u[uind:(uind + d.m - 1)])
		xind += d.n
		uind += d.m
	end
	return xarr, uarr
end

function linearize(gd::GroupDynamics, x::VF, u::VF, h::Float64)

	xarr, uarr = gsplit(gd, x, u)

	Aarr = VMF()
	Barr = VMF()
	for j = 1:gd.num_agents
		A, B = linearize(gd.array[j], xarr[j], uarr[j], h)
		push!(Aarr, A)
		push!(Barr, B)
	end

	dims = 2*ones(Int, gd.num_agents)
	dims[1] = 1

	return cat(dims, Aarr...), cat(dims, Barr...)
end


######################################################################
# integration
######################################################################
type ForwardEuler <: IntegrationScheme end
type SymplecticEuler <: IntegrationScheme end

function integrate(tm::TrajectoryManager, x::VF, u::VF)
	integrate(tm.int_scheme, tm.dynamics, x, u, tm.h)
end

function integrate(::ForwardEuler, d::Dynamics, x::VF, u::VF, h::Float64)
	forward_euler(d, x, u, h)
end

function integrate(::SymplecticEuler, d::Dynamics, x::VF, u::VF, h::Float64)
	symplectic_euler(d, x, u, h)
end

function integrate(tm::TrajectoryManager, ud::VVF)
	xd = Array(Vector{Float64}, tm.N+1)

	xd[1] = deepcopy(tm.x0)
	for i = 1:tm.N
		xd[i+1] = integrate(tm.int_scheme, tm.dynamics, x, ud[i], tm.h)
	end

	return xd
end


######################################################################
# forward_euler
######################################################################
function forward_euler(tm::TrajectoryManager, x::VF, u::VF)
	forward_euler(tm.dynamics, x, u, tm.h)
end

function forward_euler(ld::LinearDynamics, x::VF, u::VF, h::Float64)
	return ld.A*x + ld.B*u
end
function forward_euler(gd::GroupDynamics, x::VF, u::VF, h::Float64)
	n_ind = m_ind = 1
	xp = VVF()
	xarr, uarr = gsplit(gd, x, u)
	for j = 1:gd.num_agents
		push!(xp, forward_euler(gd.array[j], xarr[j], uarr[j], h))
	end
	return vcat(xp...)
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


######################################################################
# symplectic_euler
######################################################################
function symplectic_euler(d::Dynamics, x::VF, u::VF, h::Float64)
	#A,B = linearize(d, x, u, h)
	A = zeros(4,4)
	A[1,3] = 1
	A[2,4] = 1
	A[3,1] = -1
	A[4,2] = -1
	B = zeros(4,2)
	B[3,1] = 1
	B[4,2] = 1

	#n2 = round(Int, d.n/2)
	n2 = 2

	A11 = A[1:n2, 1:n2]
	A12 = A[1:n2, n2+1:d.n]
	A21 = A[n2+1:d.n, 1:n2]
	A22 = A[n2+1:d.n, n2+1:d.n]

	B1 = B[1:n2,:]
	B2 = B[n2+1:d.n,:]

	# block matrices
	Ad22 = inv(eye(n2)-h*A22)
	Ad11 = eye(n2) + h*A11 + h*h*A12*Ad22*A21
	Ad12 = h*A12*Ad22
	Ad21 = h*Ad22*A21

	Bd1 = h*h*A12*Ad22*B2 + h*B1
	Bd2 = h*Ad22*B2

	Ase = [Ad11 Ad12; Ad21 Ad22]
	Bse = [Bd1; Bd2]

	xp = Ase*x + Bse*u
	return xp
end


# split one vector per state up
function vvf2vvvf(xd::VVF, ud::VVF, vtm::Vector{TrajectoryManager})
	N = length(xd) - 1
	num_agents = length(vtm)

	x_start = u_start = 1
	xds = VVVF(num_agents)
	uds = VVVF(num_agents)
	for j = 1:num_agents
		x_end = x_start + vtm[j].dynamics.n - 1
		u_end = u_start + vtm[j].dynamics.m - 1
		xds[j] = VVF()
		uds[j] = VVF()
		for n = 1:N
			push!(xds[j], xd[n][x_start:x_end])
			push!(uds[j], ud[n][u_start:u_end])
		end
		push!(xds[j], xd[N+1][x_start:x_end])
		x_start = x_end + 1
		u_start = u_end + 1
	end
	return xds, uds
end

function vvf2vvvf(xd::VVF, vtm::Vector{TrajectoryManager})
	N = length(xd) - 1
	num_agents = length(vtm)

	x_start = 1
	xds = VVVF(num_agents)
	for j = 1:num_agents
		x_end = x_start + vtm[j].dynamics.n - 1
		xds[j] = VVF()
		for n = 1:N+1
			push!(xds[j], xd[n][x_start:x_end])
		end
		x_start = x_end + 1
	end
	return xds
end

function vvf2vvvf(xd::VVF, gd::GroupDynamics)
	N = length(xd) - 1

	x_start = 1
	xds = VVVF(gd.num_agents)
	for j = 1:gd.num_agents
		x_end = x_start + gd.array[j].n - 1
		xds[j] = VVF()
		for n = 1:(N+1)
			push!(xds[j], xd[n][x_start:x_end])
		end
		x_start = x_end + 1
	end
	return xds
end
