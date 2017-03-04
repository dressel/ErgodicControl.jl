######################################################################
# ergodicity.jl
# handles stuff needed for ergodicity
######################################################################


export ErgodicManagerSE2

"""
`ErgodicManagerSE2(L::Float64, K::Int, bins::Int)`

`ErgodicManagerSE2(example_name::String; K::Int=5, bins::Int=100)`

Valid `example_name` entries are:

* "single gaussian"
* "double gaussian"
"""
type ErgodicManagerSE2 <: ErgodicManager
	domain::Domain
	M::Int
	N::Int
	P::Int
	phi::Array{Float64,3}				# spatial distribution
	phik::Array{Complex{Float64},3}		# spatial Fourier coefficients
	Lambda::Array{Float64,3}			# constants 

	function ErgodicManagerSE2(L::Float64, K::Int, bins::Int)
		em = new()
		em.M = K
		em.N = K
		em.P = K
		em.domain = Domain([0,0,-pi/2], [L,L,pi/2], [bins,bins,bins])
		em.phik = zeros(em.M+1, em.N+1, em.P+1)
		em.Lambda = zeros(em.M+1, em.N+1, em.P+1)
		em.phi = ones(bins,bins,bins) / (bins * bins * bins)

		Lambda!(em)
		#decompose!(em)		# TODO: I should do this I think
		return em
	end

	function ErgodicManagerSE2(L::Float64, K::Int, d::Array{Float64,3})
		bins = size(d,1)	# assuming all sides have same number cells
		em = ErgodicManagerSE2(L, K, bins)
		em.phi = d

		Lambda!(em)
		decompose!(em, d)
	end

	function ErgodicManagerSE2(example_name::String; K::Int=5, bins::Int=100)
		L = 1.0
		println("h2")
		em = ErgodicManagerSE2(L, K, bins)
		println("h3")

		if example_name == "single gaussian"
			mu = [L/2.0, L/2.0, 0]
			Sigma = 0.03 * eye(3)
			println("pre fi")
			phi!(em, mu, Sigma)
			println("post fi")
			decompose!(em)
			println("post dec")
		elseif example_name == "double gaussian"
			# Create Gaussian distribution and its coefficients
			mu1 = [0.3, 0.7]
			Sigma1 = 0.025* eye(2)
			mu2 = [0.7, 0.3]
			Sigma2 = 0.025* eye(2)
			phi!(em, mu1, Sigma1, mu2, Sigma2)
			decompose!(em)
		else
			error("example name not recognized")
		end
		return em
	end
end

# fills the matrix Lambda_{m,n,p}
function Lambda!(em::ErgodicManagerSE2)
	for m = 0:em.M
		for n = 0:em.N
			for p = 0:em.P
				den_sqrt = (1.0 + m*m + n*n + p*p)
				em.Lambda[m+1, n+1, p+1] = 1.0 / (den_sqrt * den_sqrt)
			end
		end
	end
end



######################################################################
# Setting and computing phi
######################################################################
function phi!(em::ErgodicManagerSE2, dm::VF, ds::MF)
	# first, generate d
	#d = zeros(em.bins, em.bins, em.bins)
	d = zeros(x_cells(em), y_cells(em), z_cells(em))
	d_sum = 0.0
	for xi = 1:x_cells(em)
		x = x_min(em) + (xi-0.5)*x_size(em)
		println("xi = ", xi)
		for yi = 1:y_cells(em)
			y = y_min(em) + (yi-0.5)*y_size(em)
			#for zi = 1:em.domain.cells[3]
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

#function phi!(em::ErgodicManagerSE2, dm1::Vector{Float64}, ds1::Matrix{Float64}, dm2::Vector{Float64}, ds2::Matrix{Float64})
#	d = zeros(em.bins, em.bins)
#	d_sum = 0.0
#	for xi = 1:em.bins
#		x = (xi-0.5)*em.cell_size
#		for yi = 1:em.bins
#			y = (yi-0.5)*em.cell_size
#			d[xi,yi] = .5*my_pdf((x,y), dm1, ds1)+.5*my_pdf((x,y), dm2, ds2)
#			d_sum += d[xi,yi]
#		end
#	end
#	normalize!(d, em.domain.cell_size)
#	em.phi = d
#end

#function phi!(em::ErgodicManagerSE2, means::VV_F, covs::VM_F, weights::Vector{Float64})
#	half_size = em.cell_size / 2.0
#	d = zeros(em.bins, em.bins)
#	d_sum = 0.0
#	num_gaussians = length(means)
#	for xi = 1:em.bins
#		x = (xi-1)*em.cell_size + half_size
#		for yi = 1:em.bins
#			y = (yi-1)*em.cell_size + half_size
#			for gi = 1:num_gaussians
#				d[xi,yi] += weights[gi] * my_pdf((x,y), means[gi], covs[gi])
#			end
#			d_sum += d[xi,yi]
#		end
#	end
#	normalize!(d)
#	em.phi = d
#end


######################################################################
# Computing Fourier coefficients
######################################################################
# update the Fourier coefficients based on some distribution
# Here, I assume it is discrete, but I should change this...
# TODO: maybe I should do some bounds checking?
function decompose!(em::ErgodicManagerSE2, d::Array{Float64,3})
	for m = 0:em.M
		println("m = ", m)
		for n = 0:em.N
			for p = 0:em.P
				em.phik[m+1,n+1,p+1] = phi_mnp(em, m, n, p, d)
			end
		end
	end
	em.phi = d
end


# iterate over the state space
function phi_mnp(em::ErgodicManagerSE2, m::Int, n::Int, p::Int, d::Array{Float64,3})
	val = 0.0im
	#for xi = 1:em.domain.cells[1]
	for xi = 1:x_cells(em)
		x = x_min(em) + (xi-0.5) * x_size(em)
		#for yi = 1:em.domain.cells[2]
		for yi = 1:y_cells(em)
			#y = (yi-0.5) * em.domain.cell_lengths[2]
			y = y_min(em) + (yi-0.5) * y_size(em)
			#for zi = 1:em.domain.cells[3]
			for zi = 1:z_cells(em)
				z = z_min(em) + (zi-0.5) * z_size(em)
				val += d[xi,yi,zi] * F_mnp(m,n,p,x,y,z) *em.domain.cell_size
			end
		end
	end
	return val
end

function F_mnp(m::Int, n::Int, p::Int, x::Float64, y::Float64, z::Float64)
	# compute psi, r
	r = sqrt(x*x + y*y)
	psi = atan2(y,x)	# atan(y/x)
	i = float(im)
	return i^(n-m) * exp(i*(m*psi + (n-m)z)) * besselj(m-n, p*r)
end



# reconstructs from Fourier coefficients in ck
# TODO: not nearly complete
function reconstruct(em::ErgodicManagerSE2, ck::Matrix{Float64})
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
