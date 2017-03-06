######################################################################
# gaussian.jl
#
# Creates gaussian distributions over domains
######################################################################
export gaussian

gaussian(domain::Domain, dm::VF, ds::MF) = gaussian(domain, [dm], [ds])


# TODO: check that provided dimensions are correct
function gaussian(domain::Domain, means::VVF, covs::VMF, weights::VF=ones(length(means)))
	n = domain.num_dims
	if n == 2
		return gaussian2D(domain, means, covs)
	elseif n == 3
		return gaussian3D(domain, means, covs)
	end
end

function gaussian2D(domain::Domain, means::VVF, covs::VMF, weights::VF=ones(length(means)))
	d = zeros(x_cells(domain), y_cells(domain))
	d_sum = 0.0
	num_gaussians = length(means)
	for xi = 1:x_cells(domain)
		x = x_min(domain) + (xi-0.5)*x_size(domain)
		for yi = 1:y_cells(domain)
			y = y_min(domain) + (yi-0.5)*y_size(domain)
			for gi = 1:num_gaussians
				d[xi,yi] += weights[gi] * my_pdf((x,y), means[gi], covs[gi])
			end
			d_sum += d[xi,yi]
		end
	end
	normalize!(d, domain.cell_size)
	return d
end

function gaussian3D(domain::Domain, means::VVF, covs::VMF, weights::VF=ones(length(means)))
	d = zeros(x_cells(domain), y_cells(domain), z_cells(domain))
	d_sum = 0.0
	num_gaussians = length(means)
	for xi = 1:x_cells(domain)
		x = x_min(domain) + (xi-0.5)*x_size(domain)
		for yi = 1:y_cells(domain)
			y = y_min(domain) + (yi-0.5)*y_size(domain)
			for zi = 1:z_cells(domain)
				z = z_min(domain) + (zi-0.5)*z_size(domain)
				for gi = 1:num_gaussians
					d[xi,yi,zi] += weights[gi] * my_pdf([x,y,z], means[gi], covs[gi])
				end
				d_sum += d[xi,yi,zi]
			end
		end
	end
	normalize!(d, domain.cell_size)
	return d
end



######################################################################
# Old and deprecated code
#
# I am keeping it around because it might be slightly faster
#  in the event that there is only one gaussian (shouldnt matter
#   otherwise)
######################################################################

#function gaussian(domain::Domain, dm::VF, ds::MF)
#	n = length(dm)
#	@assert n == domain.num_dims
#	# TODO: check that dm and ds have correct dimensions
#
#	if n == 2
#		return gaussian2D(domain, dm, ds)
#	elseif n == 3
#		return gaussian3D(domain, dm, ds)
#	end
#end
#
#function gaussian2D(domain::Domain, dm::VF, ds::MF)
#	d = zeros(x_cells(domain), y_cells(domain))
#	d_sum = 0.0
#	for xi = 1:x_cells(domain)
#		x = x_min(domain) + (xi-0.5)*x_size(domain)
#		for yi = 1:y_cells(domain)
#			y = y_min(domain) + (yi-0.5)*y_size(domain)
#			d[xi,yi] = my_pdf((x,y), dm, ds)
#			d_sum += d[xi,yi]
#		end
#	end
#	normalize!(d, domain.cell_size)
#	return d
#end
#
#function gaussian3D(domain::Domain, dm::VF, ds::MF)
#	dims = (domain.cells...)
#	d = zeros(dims)
#	d_sum = 0.0
#	for xi = 1:x_cells(domain)
#		x = x_min(domain) + (xi-0.5)*x_size(domain)
#		for yi = 1:y_cells(domain)
#			y = y_min(domain) + (yi-0.5)*y_size(domain)
#			for zi = 1:z_cells(domain)
#				z = z_min(domain) + (zi-0.5)*z_size(domain)
#				d[xi,yi,zi] = my_pdf([x,y,z], dm, ds)
#				d_sum += d[xi,yi,zi]
#			end
#		end
#	end
#	normalize!(d, domain.cell_size)
#	return d
#end


# TODO: decide what to do with this. It is OLD
#  perhaps I should just delete it
#function phi!(em::ErgodicManagerR2, dm1::VF, ds1::MF, dm2::VF, ds2::MF)
#	d = zeros(em.domain.cells[1], em.domain.cells[2])
#	d_sum = 0.0
#	for xi = 1:em.domain.cells[1]
#		x = (xi-0.5)*em.domain.cell_lengths[1]
#		for yi = 1:em.domain.cells[2]
#			y = (yi-0.5)*em.domain.cell_lengths[2]
#			d[xi,yi] = .5*my_pdf((x,y), dm1, ds1)+.5*my_pdf((x,y), dm2, ds2)
#			d_sum += d[xi,yi]
#		end
#	end
#	normalize!(d, em.domain.cell_size)
#	em.phi = d
#end

#  dims = (domain.cells...)
#  d = zeros(d)
