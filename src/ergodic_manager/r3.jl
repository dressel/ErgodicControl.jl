######################################################################
# ergodic_manager/r3.jl
#
# handles stuff needed for ergodicity
######################################################################


export ErgodicManagerR3

"""
`ErgodicManagerR3(d::Domain, phi::Matrix{Float64}, K::Int)`
"""
type ErgodicManagerR3 <: ErgodicManager
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
	kpizl::Matrix{Float64}

	function ErgodicManagerR3(d::Domain, phi::Array{Float64,3}, K::Int=5)
		em = new()
		em.domain = deepcopy(d)
		em.K = K
		em.hk = zeros(K+1,K+1,K+1)
		em.phi = deepcopy(phi)
		em.phik = zeros(K+1,K+1,K+1)
		em.Lambda = zeros(K+1,K+1,K+1)
		em.kpixl = zeros(K+1, d.cells[1])
		em.kpiyl = zeros(K+1, d.cells[2])
		em.kpizl = zeros(K+1, d.cells[3])

		Lambda!(em)
		kpixl!(em)
		hk!(em)
		decompose!(em)
		return em
	end
end

# fills each entry Lambda[k1,k2] in the Lambda matrix
function Lambda!(em::ErgodicManagerR3)
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
function kpixl!(em::ErgodicManagerR3)
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

	Lz = em.domain.lengths[3]
	zmin = z_min(em)
	for zi = 1:z_cells(em)
		z = zmin + (zi-0.5)*z_size(em)
		for k = 0:em.K
			em.kpizl[k+1,zi] = cos(k*pi*(z-zmin) / Lz)
		end
	end
end

# generates the hk coefficients for the ergodic manager
# these coefficients only need to be computed once
function hk!(em::ErgodicManagerR3)
	for k1 = 0:em.K
		for k2 = 0:em.K
			for k3 = 0:em.K
				em.hk[k1+1,k2+1,k3+1] = hk_ij(em, k1, k2, k3)
			end
		end
	end
end

# computes the coefficients for a specific value of k1 and k2
# called by hk!
# TODO: this is wrecked by domain
function hk_ij(em::ErgodicManagerR3, k1::Int, k2::Int, k3::Int)
	val = 0.0
	for xi = 1:x_cells(em)
		cx = em.kpixl[k1+1,xi]
		cx2 = cx * cx
		for yi = 1:y_cells(em)
			cy = em.kpiyl[k2+1,yi]
			cy2 = cy * cy
			for zi = 1:z_cells(em)
				cz = em.kpizl[k3+1,zi]
				cz2 = cz * cz
				val += cx2 * cy2 * cz2 * em.domain.cell_size
			end
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
function decompose!(em::ErgodicManagerR3, d::Array{Float64,3})
	for k1 = 0:em.K
		for k2 = 0:em.K
			for k3 = 0:em.K
				em.phik[k1+1,k2+1,k3+1] = phi_ij(em, k1, k2, k3, d)
			end
		end
	end
	em.phi = d
end


# iterate over the state space
function phi_ij(em::ErgodicManagerR3, k1::Int, k2::Int, k3::Int, d::Array{Float64,3})
	val = 0.0
	for xi = 1:x_cells(em)
		cx = em.kpixl[k1+1,xi]
		for yi = 1:y_cells(em)
			cy = em.kpiyl[k2+1,yi]
			for zi = 1:z_cells(em)
				cz = em.kpizl[k3+1,zi]
				val += d[xi,yi,zi] * cx * cy * cz * em.domain.cell_size
			end
		end
	end
	return val / em.hk[k1+1,k2+1,k3+1]
end


function decompose(em::ErgodicManagerR3, traj::VVF)
	N = length(traj)-1

	# Array to hold trajectory's Fourier coefficients
	ck = zeros(em.K+1, em.K+1, em.K+1)

	# length of each dimension
	Lx = em.domain.lengths[1]
	Ly = em.domain.lengths[2]
	Lz = em.domain.lengths[3]

	# minimum values in each dimension
	xmin = x_min(em)
	ymin = y_min(em)
	zmin = z_min(em)

	for k1 = 0:em.K
		kpiL1 = k1 * pi / Lx
		for k2 = 0:em.K
			kpiL2 = k2 * pi / Ly
			for k3 = 0:em.K
				kpiL3 = k3 * pi / Lz
				hk = em.hk[k1+1, k2+1, k3+1]
				fk_sum = 0.0
				# now loop over time
				for n = 0:N
					xn = traj[n + 1]
					c1 = cos(kpiL1 * (xn[1] - xmin))
					c2 = cos(kpiL2 * (xn[2] - ymin))
					c3 = cos(kpiL3 * (xn[3] - zmin))
					fk_sum += c1*c2*c3
				end
				ck[k1+1, k2+1, k3+1] = fk_sum / (hk * (N+1))
			end
		end
	end
	return ck
end
