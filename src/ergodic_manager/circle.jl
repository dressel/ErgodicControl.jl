######################################################################
# circle.jl
######################################################################

export circle

"""
`circle(d::Domain, xc::VF, r::Float64)`

Returns a distribution that evenly weights points in a circle centered at `xc` with radius `r`. The distribution is normalized (so it integrates to 1) before being returned to the user.
"""
function circle(d::Domain, xc::Vector{Float64}, radius::Float64)
	phi = zeros(x_cells(d), y_cells(d))
	R2 = radius * radius

	for xi = 1:x_cells(d)
		x = x_min(d) + (xi-0.5)*x_size(d)
		for yi = 1:y_cells(d)
			y = y_min(d) + (yi-0.5)*y_size(d)
			dx = xc[1] - x
			dy = xc[2] - y
			if (dx*dx + dy*dy) < R2
				phi[xi,yi] = 1.0
			end
		end
	end
	normalize!(phi, d.cell_size)
	return phi
end
