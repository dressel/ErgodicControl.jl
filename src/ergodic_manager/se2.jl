######################################################################
# se2.jl
#
######################################################################

const enum = enumerate


export ErgodicManagerSE2

"""
`ErgodicManagerSE2(L::Float64, K::Int, bins::Int)`

`ErgodicManagerSE2(example_name::String; K::Int=5, bins::Int=100)`

Valid `example_name` entries are:

* "single gaussian"
* "double gaussian"
"""
mutable struct ErgodicManagerSE2 <: ErgodicManager
    domain::Domain

    M::UnitRange{Int}               # Fourier coefficient indices
    N::UnitRange{Int}
    P::Int

    phi::Array{Float64,3}				# spatial distribution
    phik::Array{Complex{Float64},3}		# spatial Fourier coefficients

    Lambda::Array{Float64,3}			# constants 
    cache::Array{Complex{Float64},6}    # stores values of F_mnp(x,y,z)


    function ErgodicManagerSE2(d::Domain, phi::Array{Float64,3}, K::Int=5)
        em = new()
        em.domain = deepcopy(d)

        em.M = 0:K
        em.N = 0:K
        em.P = K

        em.phi = deepcopy(phi)
        em.phik = zeros(length(em.M), length(em.N), em.P)
        em.Lambda = zeros(length(em.M), length(em.N), em.P)

        tic()
        cache!(em)
        println("SE2 cache time = ", round(toq(),3), " seconds")

        Lambda!(em)
        decompose!(em)
        return em
    end

    function ErgodicManagerSE2(d::Domain, K::Int=5)
        em = new()
        em.domain = deepcopy(d)

        em.M = 0:K
        em.N = 0:K
        em.P = K

        em.phik = zeros(em.M+1, em.N+1, em.P)
        em.Lambda = zeros(em.M+1, em.N+1, em.P)

        Lambda!(em)
        return em
    end
end


# Calls to F_mnp are expensive as they rely on the bessel function.
# Therefore, we store calls to this function in an array beforehand.
function cache!(em::ErgodicManagerSE2)
    n_M = length(em.M)
    n_N = length(em.N)
    cache = zeros(Complex{Float64}, z_cells(em), y_cells(em), x_cells(em), em.P, n_N, n_M)
    for (mi,m) in enum(em.M), (ni,n) in enum(em.N), p = 1:em.P
        for xi = 1:x_cells(em), yi = 1:y_cells(em), zi=1:z_cells(em)
            x = x_min(em) + (xi-0.5)*x_size(em)
            y = y_min(em) + (yi-0.5)*y_size(em)
            z = z_min(em) + (zi-0.5)*z_size(em)
            cache[zi,yi,xi,p,ni,mi] = F_mnp(m,n,p,x,y,z)
        end
    end
    em.cache = cache
end


# fills the matrix Lambda_{m,n,p}
function Lambda!(em::ErgodicManagerSE2)
    for (mi,m) in enum(em.M), (ni,n) in enum(em.N), p = 1:em.P
        den_sqrt = (1.0 + m*m + n*n + p*p)
        em.Lambda[mi, ni, p] = 1.0 / (den_sqrt * den_sqrt)
    end
end


######################################################################
# Setting and computing phi
######################################################################
function phi!(em::ErgodicManagerSE2, dm::VF, ds::MF)
    # first, generate d
    d = zeros(x_cells(em), y_cells(em), z_cells(em))
    d_sum = 0.0
    for xi = 1:x_cells(em)
		x = x_min(em) + (xi-0.5)*x_size(em)
		println("xi = ", xi)
		for yi = 1:y_cells(em)
			y = y_min(em) + (yi-0.5)*y_size(em)
			for zi = 1:z_cells(em)
				z = z_min(em) + (zi-0.5)*z_size(em)
				d[xi,yi,zi] = my_pdf([x,y,z], dm, ds)
				d_sum += d[xi,yi,zi]
			end
		end
	end
    normalize!(d, em.domain.cell_size)
    em.phi = d
end


######################################################################
# Computing Fourier coefficients
######################################################################
# update the Fourier coefficients based on some distribution
# Here, I assume it is discrete, but I should change this...
# TODO: maybe I should do some bounds checking?
function decompose!(em::ErgodicManagerSE2, d::Array{Float64,3})
    for (mi,m) in enum(em.M), (ni,n) in enum(em.N), p = 1:em.P
        em.phik[mi,ni,p] = phi_mnp(em, mi, ni, p, d)
    end
    em.phi = d
end


# computes phi_mnp by iterating over the state space
# uses the cache to evaluate F_mnp
function phi_mnp(em::ErgodicManagerSE2, mi::Int, ni::Int, p::Int, d::Array{Float64,3})
	val = 0.0im
	for xi = 1:x_cells(em), yi = 1:y_cells(em), zi = 1:z_cells(em)
        val += d[xi,yi,zi] * em.cache[zi, yi, xi, p, ni, mi]
	end
	return val * em.domain.cell_size
end

function F_mnp(m::Int, n::Int, p::Int, x::Float64, y::Float64, z::Float64)
    # compute psi, r
    r = sqrt(x*x + y*y)
    psi = atan2(y,x)	# atan(y/x)
    i = float(im)
    return i^(n-m) * exp(i*(m*psi + (n-m)z)) * besselj(m-n, p*r)
end


function decompose(em::ErgodicManagerSE2, traj::VVF)
    N = length(traj)-1

    ck = zeros(Complex{Float64}, length(em.M), length(em.N), em.P)

    for (mi,m) in enum(em.M), (ni,n) in enum(em.N), p = 1:em.P
        fk_sum = 0.0im
        # now loop over time
        for i = 0:N-1
            xi = traj[i + 1]
            fk_sum += F_mnp(m, n, p, xi[1], xi[2], xi[3])
        end
        # TODO: check that this is right
        ck[mi, ni, p] = fk_sum / N
    end
    return ck
end


# reconstructs from Fourier coefficients in ck
# This comes from last equation in chapter 10.3 of 
#  Engineering Applications of Noncommutative Harmonic Analysis
#
# Returns an array of complex numbers.
# Even if you expect an array of reals, there will be small complex values
# call real(output) if you want just the real parts
function reconstruct(em::ErgodicManagerSE2, ck::Array{Complex{Float64},3})

    vals = zeros(Complex{Float64}, x_cells(em), y_cells(em), z_cells(em))

    # this order is much better
    for (mi,m) in enum(em.M), (ni,n) in enum(em.N), p = 1:em.P
        c = ck[mi, ni, p]
        for xi = 1:x_cells(em), yi = 1:y_cells(em), zi=1:z_cells(em)
            f = conj(em.cache[zi, yi, xi, p, ni, mi])
            vals[xi,yi,zi] += c * f * p
        end
	end
	return vals
end
