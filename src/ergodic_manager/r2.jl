######################################################################
# ergodic_manager/r2.jl
#
# handles stuff needed for ergodicity
######################################################################


export ErgodicManagerR2

# kpixl is a (K+1 x bins) matrix, each entry storing cos(k*pi*x / L)
#  this assumes some discretization
"""
`ErgodicManager(d::Domain, phi::Matrix{Float64}, K::Int)`

`ErgodicManager(example_name::String; K::Int=5, bins::Int=100)`

Valid `example_name` entries are:

* "single gaussian"
* "double gaussian"
"""
type ErgodicManagerR2 <: ErgodicManager
	domain::Domain				# spatial domain
	K::Int						# number of Fourier coefficients
	phi::Matrix{Float64}		# spatial distribution
	phik::Matrix{Float64}		# distribution's Fourier coefficients

	# constant regardless of phi (depend on k1,k2)
	Lambda::Matrix{Float64}
	hk::Matrix{Float64}

	# to speed up computation
	kpixl::Matrix{Float64}
	kpiyl::Matrix{Float64}

	function ErgodicManagerR2(d::Domain, phi::MF, K::Int=5)
		em = new()
		em.domain = deepcopy(d)
		em.K = K
		em.hk = zeros(K+1,K+1)
		em.phi = deepcopy(phi)
		em.phik = zeros(K+1,K+1)
		em.Lambda = zeros(K+1,K+1)
		em.kpixl = zeros(K+1, d.cells[1])
		em.kpiyl = zeros(K+1, d.cells[2])

		Lambda!(em)
		kpixl!(em)
		hk!(em)
		decompose!(em)
		return em
	end

	function ErgodicManagerR2(d::Domain, K::Int=5)
		em = new()
		em.domain = deepcopy(d)
		em.K = K
		em.hk = zeros(K+1,K+1)
		#em.phi = deepcopy(phi)
		em.phik = zeros(K+1,K+1)
		em.Lambda = zeros(K+1,K+1)
		em.kpixl = zeros(K+1, d.cells[1])
		em.kpiyl = zeros(K+1, d.cells[2])

		Lambda!(em)
		kpixl!(em)
		hk!(em)
		#decompose!(em)
		return em
	end
end


# fills each entry Lambda[k1,k2] in the Lambda matrix
function Lambda!(em::ErgodicManagerR2)
    for k1 = 0:em.K, k2 = 0:em.K
        den = (1.0 + k1*k1 + k2*k2) ^ 1.5
        em.Lambda[k1+1, k2+1] = 1.0 / den
    end
end


# TODO: check if the row/col ordering of this is julia-efficient
# TODO: also, I could probably speed up by iterating over k first?
function kpixl!(em::ErgodicManagerR2)
	Lx = em.domain.lengths[1]
	xmin = x_min(em)
	for xi = 1:x_cells(em)
		# TODO: seems like I just cancel xmin out. Ignore it maybe?
		x = xmin + (xi-0.5)*x_size(em)
		for k = 0:em.K
			em.kpixl[k+1,xi] = cos(k*pi*(x-xmin) / Lx)
		end
	end

	Ly = em.domain.lengths[2]
	ymin = y_min(em)
	for yi = 1:y_cells(em)
		y = ymin + (yi-0.5)*y_size(em)
		for k = 0:em.K
			em.kpiyl[k+1,yi] = cos(k*pi*(y-ymin) / Ly)
		end
	end
end

# generates the hk coefficients for the ergodic manager
# these coefficients only need to be computed once
function hk!(em::ErgodicManagerR2)
    for k1 = 0:em.K, k2 = 0:em.K
        em.hk[k1+1,k2+1] = hk_ij(em, k1, k2)
    end
end

# computes the coefficients for a specific value of k1 and k2
# called by hk!
function hk_ij(em::ErgodicManagerR2, k1::Int, k2::Int)
	val = 0.0
	for xi = 1:x_cells(em)
		cx = em.kpixl[k1+1,xi]
		cx2 = cx * cx
		for yi = 1:y_cells(em)
			cy = em.kpiyl[k2+1,yi]
			val += cx2 * cy * cy * em.domain.cell_size
		end
	end

	return sqrt(val)
end



######################################################################
# Computing Fourier coefficients
######################################################################
# update the Fourier coefficients based on some distribution
# Here, I assume it is discrete, but I should change this...
# TODO: maybe I should do some bounds checking?
function decompose!(em::ErgodicManagerR2, d::Matrix{Float64})
    for k1 = 0:em.K, k2 = 0:em.K
        em.phik[k1+1,k2+1] = phi_ij(em, k1, k2, d)
    end
    em.phi = d
end


# iterate over the state space
function phi_ij(em::ErgodicManagerR2, k1::Int, k2::Int, d::Matrix{Float64})
	val = 0.0
	for xi = 1:x_cells(em)
		cx = em.kpixl[k1+1,xi]
		for yi = 1:y_cells(em)
			cy = em.kpiyl[k2+1,yi]
			val += d[xi,yi] * cx * cy * em.domain.cell_size
		end
	end
	return val / em.hk[k1+1,k2+1]
end

"""
`decompose(em, traj::VVF)`

Decomposes a set of positions into a set of `ck` Fourier coefficients.

"""
function decompose(em::ErgodicManagerR2, traj::VVF)

	# trajectory is really of length N+1
	N = length(traj)-1

	# create matrix to hold trajectory's Fourier coefficients
	ck = zeros(em.K+1, em.K+1)

	# lengths of each dimension
	Lx = em.domain.lengths[1]
	Ly = em.domain.lengths[2]

	# minimum values in each dimension
	xmin = x_min(em)
	ymin = y_min(em)

	for k1 = 0:em.K
		kpiL1 = k1 * pi / Lx
		for k2 = 0:em.K
			kpiL2 = k2 * pi / Ly
			hk = em.hk[k1+1, k2+1]
			fk_sum = 0.0
			# now loop over time
			for n = 0:N
				xn = traj[n + 1]
				c1 = cos(kpiL1 * (xn[1]-xmin))
				c2 = cos(kpiL2 * (xn[2]-ymin))
				fk_sum += c1*c2
			end
			ck[k1+1, k2+1] = fk_sum / (hk * (N+1))
		end
	end
	return ck
end


# reconstructs from Fourier coefficients in ck
function reconstruct(em::ErgodicManagerR2, ck::Matrix{Float64})
	# iterate over all bins
	v = em.domain.cell_size
	vals = zeros(x_cells(em), y_cells(em))

	for xi = 1:x_cells(em)
		x = x_min(em) + (xi-0.5)*x_size(em)
		for yi = 1:y_cells(em)
			y = y_min(em) + (yi-0.5)*y_size(em)
			for k1 = 0:em.K
				cx = em.kpixl[k1+1,xi]
				for k2 = 0:em.K
					cy = em.kpixl[k2+1,yi]
					vals[xi,yi] += ck[k1+1,k2+1]*cx*cy/em.hk[k1+1,k2+1]
				end
			end
		end
	end
	return vals
end
