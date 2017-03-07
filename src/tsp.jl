######################################################################
# tsp.jl
#
# Some code for the traveling salesman problem
######################################################################

# uses nearest neighbors to heuristically solve tsp
# tsp is traveling salesman problem
function tsp_nn!(xd::VVF, points::VVF)
	xc = xd[1]
	next_p = xd[1]
	n = 1
	while length(points) > 0
		best_d = Inf
		best_ind = 1
		for (p_ind,p) in enumerate(points)
			r = xc - p 
			d = dot(r, r)
			if d < best_d
				best_d = d
				best_ind = p_ind
			end
		end
		xd[n+1] = points[best_ind]
		deleteat!(points, best_ind)
		n += 1
	end
end


# simply takes the points as they are
# this creates a shitty trajectory
function tsp_rand!(xd::VVF, points::VVF)
	n = 1
	for p in points
		xd[n+1] = p
		n += 1
	end
end
