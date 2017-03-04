######################################################################
# gaussian.jl
#
# Creates gaussian distributions over domains
######################################################################
export gaussian

# TODO: handle the three dimensional case as well
function gaussian(domain::Domain, dm::VF, ds::MF)
	d = zeros(x_cells(domain), y_cells(domain))
	d_sum = 0.0
	for xi = 1:x_cells(domain)
		x = x_min(domain) + (xi-0.5)*x_size(domain)
		for yi = 1:y_cells(domain)
			y = y_min(domain) + (yi-0.5)*y_size(domain)
			d[xi,yi] = my_pdf((x,y), dm, ds)
			d_sum += d[xi,yi]
		end
	end
	normalize!(d, domain.cell_size)
	return d
end



#function phi!(em::ErgodicManagerR2, dm::VF, ds::MF)
#	# first, generate d
#	d = zeros(x_cells(em), y_cells(em))
#	d_sum = 0.0
#	for xi = 1:x_cells(em)
#		x = x_min(em) + (xi-0.5)*x_size(em)
#		for yi = 1:y_cells(em)
#			y = y_min(em) + (yi-0.5)*y_size(em)
#			d[xi,yi] = my_pdf((x,y), dm, ds)
#			d_sum += d[xi,yi]
#		end
#	end
#	normalize!(d, em.domain.cell_size)
#	em.phi = d
#end

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

#function phi!(em::ErgodicManagerR2, means::VVF, covs::VMF, weights::VF)
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
