######################################################################
# ergodicity.jl
# handles stuff needed for ergodicity
######################################################################


# kpixl is a (K+1 x bins) matrix, each entry storing cos(k*pi*x / L)
#  this assumes some discretization
"""
`ErgodicManager(L::Float64, K::Int, bins::Int)`

`ErgodicManager(example_name::String; K::Int=5, bins::Int=100)`

Valid `example_name` entries are:

* "single gaussian"
* "double gaussian"
"""
type ErgodicManager
	K::Int
	bins::Int				# number of bins per side
	L::Float64
	cell_size::Float64
	hk::Matrix{Float64}
	phi::Matrix{Float64}
	phik::Matrix{Float64}
	Lambdak::Matrix{Float64}
	kpixl::Matrix{Float64}

	function ErgodicManager(L::Float64, K::Int, bins::Int)
		hk   = zeros(K+1,K+1)
		phik = zeros(K+1,K+1)
		Lambdak = zeros(K+1,K+1)
		kpixl = zeros(K+1, bins)
		cell_size = L / bins
		phi = ones(bins,bins) / (bins * bins)
		em = new(K, bins, L, cell_size, hk, phi, phik, Lambdak, kpixl)

		Lambdak!(em)
		kpixl!(em)
		hk!(em)
		return em
	end

	function ErgodicManager(L::Float64, K::Int, d::Matrix{Float64})
		hk   = zeros(K+1,K+1)
		phik = zeros(K+1,K+1)
		Lambdak = zeros(K+1,K+1)
		kpixl = zeros(K+1, bins)
		cell_size = L / bins
		#phi = ones(bins,bins) / (bins * bins)
		em = new(K, bins, L, cell_size, hk, d, phik, Lambdak, kpixl)

		Lambdak!(em)
		kpixl!(em)
		hk!(em)
		phik!(em, d)
	end

	function ErgodicManager(example_name::String; K::Int=5, bins::Int=100)
		L = 1.0
		em = ErgodicManager(L, K, bins)

		if example_name == "single gaussian"
			mu = [L/2.0, L/2.0]
			Sigma = 0.03 * eye(2)
			phik!(em, mu, Sigma)
		elseif example_name == "double gaussian"
			# Create Gaussian distribution and its coefficients
			mu1 = [0.3, 0.7]
			Sigma1 = 0.025* eye(2)
			mu2 = [0.7, 0.3]
			Sigma2 = 0.025* eye(2)
			phik!(em, mu1, Sigma1, mu2, Sigma2)
		else
			error("example name not recognized")
		end
		return em
	end
end

# fills the Alpha_{k1,k2}
function Lambdak!(em::ErgodicManager)
	for k1 = 0:em.K
		for k2 = 0:em.K
			den = (1.0 + k1*k1 + k2*k2) ^ 1.5
			em.Lambdak[k1+1, k2+1] = 1.0 / den
		end
	end
end

# TODO: check if the row/col ordering of this is julia-efficient
function kpixl!(em::ErgodicManager)
	half_size = em.cell_size / 2.0
	for xi = 1:em.bins
		x = (xi-1)*em.cell_size + half_size
		for k = 0:em.K
			em.kpixl[k+1,xi] = cos(k*pi*x / em.L)
		end
	end
end

# generates the hk coefficients for the ergodic manager
# these coefficients only need to be computed once
function hk!(em::ErgodicManager)
	for k1 = 0:em.K
		for k2 = 0:em.K
			em.hk[k1+1,k2+1] = hk_ij(em, k1, k2)
		end
	end
end

# computes the coefficients for a specific value of k1 and k2
# called by hk!
function hk_ij(em::ErgodicManager, k1::Int, k2::Int)
	cs2 = em.cell_size * em.cell_size
	val = 0.0
	for xi = 1:em.bins
		cx = em.kpixl[k1+1,xi]
		cx2 = cx * cx
		for yi = 1:em.bins
			cy = em.kpixl[k2+1,yi]
			val += cx2 * cy * cy * cs2
		end
	end

	return sqrt(val)
end


######################################################################
# Setting and computing phi
######################################################################
function phi!(em::ErgodicManager, dm::Vector{Float64}, ds::Matrix{Float64})
	# first, generate d
	half_size = em.cell_size / 2.0
	d = zeros(em.bins, em.bins)
	d_sum = 0.0
	for xi = 1:em.bins
		x = (xi-1)*em.cell_size + half_size
		for yi = 1:em.bins
			y = (yi-1)*em.cell_size + half_size
			d[xi,yi] = my_pdf((x,y), dm, ds)
			d_sum += d[xi,yi]
		end
	end
	normalize!(d)
	em.phi = d
end

function phi!(em::ErgodicManager, dm1::Vector{Float64}, ds1::Matrix{Float64}, dm2::Vector{Float64}, ds2::Matrix{Float64})
	half_size = em.cell_size / 2.0
	d = zeros(em.bins, em.bins)
	d_sum = 0.0
	for xi = 1:em.bins
		x = (xi-1)*em.cell_size + half_size
		for yi = 1:em.bins
			y = (yi-1)*em.cell_size + half_size
			d[xi,yi] = .5*my_pdf((x,y), dm1, ds1)+.5*my_pdf((x,y), dm2, ds2)
			d_sum += d[xi,yi]
		end
	end
	normalize!(d)
	em.phi = d
end

function phi!(em::ErgodicManager, means::VV_F, covs::VM_F, weights::Vector{Float64})
	half_size = em.cell_size / 2.0
	d = zeros(em.bins, em.bins)
	d_sum = 0.0
	num_gaussians = length(means)
	for xi = 1:em.bins
		x = (xi-1)*em.cell_size + half_size
		for yi = 1:em.bins
			y = (yi-1)*em.cell_size + half_size
			for gi = 1:num_gaussians
				d[xi,yi] += weights[gi] * my_pdf((x,y), means[gi], covs[gi])
			end
			d_sum += d[xi,yi]
		end
	end
	normalize!(d)
	em.phi = d
end


######################################################################
# Computing Fourier coefficients
######################################################################
# update the Fourier coefficients based on some distribution
# Here, I assume it is discrete, but I should change this...
# TODO: maybe I should do some bounds checking?
function phik!(em::ErgodicManager, d::Matrix{Float64})
	for k1 = 0:em.K
		for k2 = 0:em.K
			em.phik[k1+1,k2+1] = phik_ij(em, k1, k2, d)
		end
	end
	em.phi = d
end

function phik!(em::ErgodicManager, dm::Vector{Float64}, ds::Matrix{Float64})
	phi!(em, dm, ds)
	phik!(em)
end

function phik!(em::ErgodicManager, dm1::Vector{Float64}, ds1::Matrix{Float64}, dm2::Vector{Float64}, ds2::Matrix{Float64})
	phi!(em, dm1, ds1, dm2, ds2)
	phik!(em)
end

function phik!(em::ErgodicManager, means::VV_F, covs::VM_F, weights::Vector{Float64})
	phi!(em, means, covs, weights)
	phik!(em)
end

phik!(em::ErgodicManager) = phik!(em, em.phi)


# iterate over the state space
function phik_ij(em::ErgodicManager, k1::Int, k2::Int, d::Matrix{Float64})
	val = 0.0
	cs2 = em.cell_size * em.cell_size
	for xi = 1:em.bins
		cx = em.kpixl[k1+1,xi]
		for yi = 1:em.bins
			cy = em.kpixl[k2+1,yi]
			val += d[xi,yi] * cx * cy * cs2
		end
	end
	return val / em.hk[k1+1,k2+1]
end



"""
`phi = reconstruct(em)`

Reconstructs distribution from Fourier coefficients `em.phik`.

`phi = reconstruct(em, ck::Matrix{Float64})`

Reconstructs distribution from Fourier coefficients `ck`.

Both options return a distribution `phi`, a `em.bins` by `em.bins` matrix.
The value at `phi[xi,yi]` gives the value at the cell index `xi,yi`.

`phi = reconstruct(em, xd::VV_F, start_idx=1)`

Reconstructs from trajectory `xd`. If `start_idx` is 1, then the right Riemann sum is used. If it is 0, the left Riemann sum is used.
"""
reconstruct(em::ErgodicManager) = reconstruct(em, em.phik)

# reconstructs from Fourier coefficients in ck
function reconstruct(em::ErgodicManager, ck::Matrix{Float64})
	# iterate over all bins
	half_size = em.cell_size / 2.0
	cs2 = em.cell_size * em.cell_size

	vals = zeros(em.bins, em.bins)

	for xi = 1:em.bins
		x = (xi-1)*em.cell_size + half_size
		for yi = 1:em.bins
			y = (yi-1)*em.cell_size + half_size
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
