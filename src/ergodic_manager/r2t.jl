######################################################################
# ergodic_manager/r2t.jl
#
# Euclidean (2d), with time-dependence
######################################################################


export ErgodicManagerR2T

# kpixl is a (K+1 x bins) matrix, each entry storing cos(k*pi*x / L)
#  this assumes some discretization
"""
`ErgodicManagerR2T(d::Domain, phi::Matrix{Float64}, K::Int)`
"""
type ErgodicManagerR2T <: ErgodicManager
	domain::Domain				# spatial domain
	K::Int						# number of Fourier coefficients
	phi::Array{Float64,3}		# spatial distribution
	phik::Array{Float64,3}		# distribution's Fourier coefficients

	# constant regardless of phi (depend on k1,k2)
	Lambda::Array{Float64,3}
	hk::Array{Float64,3}

	# to speed up computation
	kpixl::Matrix{Float64}
	kpiyl::Matrix{Float64}
	kpitl::Matrix{Float64}

	function ErgodicManagerR2T(d::Domain, phi::Array{Float64,3}, K::Int=5)
		em = new()
		em.domain = deepcopy(d)
		em.K = K
		em.hk = zeros(K+1,K+1,K+1)
		em.phi = deepcopy(phi)
		em.phik = zeros(K+1,K+1,K+1)
		em.Lambda = zeros(K+1,K+1,K+1)
		em.kpixl = zeros(K+1, d.cells[1])
		em.kpiyl = zeros(K+1, d.cells[2])
		N = size(phi, 3) - 1
		em.kpitl = zeros(K+1, N+1)

		Lambda!(em)
		kpixl!(em, N)
		hk!(em, N)
		decompose!(em)
		return em
	end


end

# fills each entry Lambda[k1,k2] in the Lambda matrix
# TODO: this is just repeated stuff, I should make this better
function Lambda!(em::ErgodicManagerR2T)
	for k1 = 0:em.K
		for k2 = 0:em.K
			for k3 = 0:em.K
				den = (1.0 + k1*k1 + k2*k2 + k3*k3) ^ 2.0
				em.Lambda[k1+1, k2+1, k3+1] = 1.0 / den
			end
		end
	end
end


# TODO: check if the row/col ordering of this is julia-efficient
# TODO: also, I could probably speed up by iterating over k first?
function kpixl!(em::ErgodicManagerR2T, N::Int)
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

	# cos(k*pi*t/T) = cos(k*pi*(n*dt) / (N*dt)) = cos(k*pi*n / N)
	for n = 0:N
		for k = 0:em.K
			em.kpitl[k+1,n+1] = cos(k*pi*n / N)
		end
	end
end

# generates the hk coefficients for the ergodic manager
# these coefficients only need to be computed once
function hk!(em::ErgodicManagerR2T, N::Int)
	for k1 = 0:em.K
		for k2 = 0:em.K
			for k3 = 0:em.K
				em.hk[k1+1,k2+1,k3+1] = hk_ij(em, k1, k2, k3, N)
			end
		end
	end
end

# computes the coefficients for a specific value of k1 and k2
# called by hk!
function hk_ij(em::ErgodicManagerR2T, k1::Int, k2::Int, k3::Int, N::Int)
	val = 0.0
	for xi = 1:x_cells(em)
		cx = em.kpixl[k1+1,xi]
		cx2 = cx * cx
		for yi = 1:y_cells(em)
			cy = em.kpiyl[k2+1,yi]
			cy2 = cy * cy
			for n = 0:N
				ct = em.kpitl[k3+1, n+1]
				ct2 = ct * ct
				val += cx2 * cy2 * ct2 * em.domain.cell_size
			end
		end
	end

	return sqrt(val / (N+1))
end



######################################################################
# Computing Fourier coefficients
######################################################################
# update the Fourier coefficients based on some distribution
# Here, I assume it is discrete, but I should change this...
# TODO: maybe I should do some bounds checking?
function decompose!(em::ErgodicManagerR2T, d::Array{Float64,3})
	N = size(d, 3) - 1
	for k1 = 0:em.K
		for k2 = 0:em.K
			for k3 = 0:em.K
				em.phik[k1+1,k2+1,k3+1] = phi_ij(em, k1, k2, k3, d,N)/(N+1)
			end
		end
	end
	em.phi = d
end


# iterate over the state space
function phi_ij(em::ErgodicManagerR2T, k1::Int, k2::Int, k3::Int, d::Array{Float64,3}, N::Int)
	val = 0.0
	for xi = 1:x_cells(em)
		cx = em.kpixl[k1+1,xi]
		for yi = 1:y_cells(em)
			cy = em.kpiyl[k2+1,yi]
			for n = 0:N
				ct = em.kpitl[k3+1,n+1]
				val += d[xi,yi,n+1] * cx * cy * ct * em.domain.cell_size
			end
		end
	end
	return val / em.hk[k1+1,k2+1,k3+1]
end


# reconstructs from Fourier coefficients in ck
#function reconstruct(em::ErgodicManagerR2, ck::Matrix{Float64})
#	# iterate over all bins
#	half_size = em.cell_size / 2.0
#	cs2 = em.cell_size * em.cell_size
#
#	vals = zeros(em.bins, em.bins)
#
#	for xi = 1:em.bins
#		x = (xi-1)*em.cell_size + half_size
#		for yi = 1:em.bins
#			y = (yi-1)*em.cell_size + half_size
#			for k1 = 0:em.K
#				cx = em.kpixl[k1+1,xi]
#				for k2 = 0:em.K
#					cy = em.kpixl[k2+1,yi]
#					vals[xi,yi] += ck[k1+1,k2+1]*cx*cy/em.hk[k1+1,k2+1]
#				end
#			end
#		end
#	end
#	return vals
#end
