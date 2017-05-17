######################################################################
# gaussian.jl
#
# Creates gaussian distributions over domains
######################################################################
export gaussian

gaussian(domain::Domain, dm::VF, ds::MF) = gaussian(domain, [dm], [ds])


"""
`gaussian(domain, means, covs, weights)`

Generates a distribution over `domain` from a weighted collection of Gaussians. The `means` argument is a vector of mean vectors; `covs` is a vector of covariance matrices; `weights` is a vector of real-valued numbers weighting each Gaussian. The weights need not add to 1, and the argument is optional; if `weights` is not provided, all Gaussians are weighted equally.

The resulting distribution is normalized.

The simple constructor `gaussian(domain, mean, cov)` allows you to construct a distribution from a single Gaussian.
"""
function gaussian{T<:Real}(domain::Domain, means::VVF, covs::VMF, weights::Vector{T}=ones(length(means)))
	# TODO: check that provided dimensions are correct
	w = float(weights)
	n = domain.num_dims
	if n == 2
		return gaussian2D(domain, means, covs, w)
	elseif n == 3
		return gaussian3D(domain, means, covs, w)
	end
end

function gaussian2D(domain::Domain, means::VVF, covs::VMF, weights::VF)
	d = zeros(x_cells(domain), y_cells(domain))
	d_sum = 0.0
	num_gaussians = length(means)

	inv_covs, dens = get_fast_tools(means, covs)

	for xi = 1:x_cells(domain)
		x = x_min(domain) + (xi-0.5)*x_size(domain)
		for yi = 1:y_cells(domain)
			y = y_min(domain) + (yi-0.5)*y_size(domain)
			for gi = 1:num_gaussians
				d[xi,yi] += weights[gi] * fast_pdf([x,y], means[gi], inv_covs[gi], dens[gi])
			end
			d_sum += d[xi,yi]
		end
	end
	normalize!(d, domain.cell_size)
	return d
end

function gaussian3D(domain::Domain, means::VVF, covs::VMF, weights::VF)
	d = zeros(x_cells(domain), y_cells(domain), z_cells(domain))
	d_sum = 0.0
	num_gaussians = length(means)

	# create inv_covs matrix, matrix of inverse covariants
	inv_covs, dens = get_fast_tools(means, covs)

	for xi = 1:x_cells(domain)
		x = x_min(domain) + (xi-0.5)*x_size(domain)
		for yi = 1:y_cells(domain)
			y = y_min(domain) + (yi-0.5)*y_size(domain)
			for zi = 1:z_cells(domain)
				z = z_min(domain) + (zi-0.5)*z_size(domain)
				for gi = 1:num_gaussians
					d[xi,yi,zi] += weights[gi] * fast_pdf([x,y,z], means[gi], inv_covs[gi], dens[gi])
				end
				d_sum += d[xi,yi,zi]
			end
		end
	end
	normalize!(d, domain.cell_size)
	return d
end

# Must compute Gaussian everywhere in domain
# To speed this up, we precompute inverse of cov matrix and denominator
# Each of these take a long time
function get_fast_tools(means::VVF, covs::VMF)
	num_gaussians = length(means)
	k = length(means[1])

	inv_covs = deepcopy(covs)
	dens = zeros(num_gaussians)
	for gi = 1:num_gaussians
		inv_covs[gi] = inv(covs[gi])
		dens[gi] = sqrt( (2.0pi)^k * det(covs[gi]) )
	end
	return inv_covs, dens
end
