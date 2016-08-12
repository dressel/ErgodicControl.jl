######################################################################
# includes stuff I don't want to need distributions for...
######################################################################


# only handles stuff in 2d
function my_pdf(x::Vector{Float64}, m::Vector{Float64}, S::Matrix{Float64})
	return my_pdf((x[1],x[2]), m, S)
end

function my_pdf(x::NTuple{2,Float64},m::Vector{Float64},S::Matrix{Float64})
	dx = x[1] - m[1]
	dy = x[2] - m[2]

	Si = inv(S)
	e_stuff = dx*Si[1,1]*dx + dx*Si[1,2]*dy + dy*Si[2,1]*dx + dy*Si[2,2]*dy
	return det(sqrtm(Si)) * exp(-0.5 * e_stuff) / (2.0 * pi)
end

function centroid(d::Matrix{Float64}, L::Float64)
	x_val = 0.0; y_val = 0.0
	n = size(d,1)
	cell_size = L / n
	d_sum = 0.0
	for x = 1:n
		for y = 1:n
			x_val += (x-.5) * d[x,y]
			y_val += (y-.5) * d[x,y]
			d_sum += d[x,y]
		end
	end
	return x_val*cell_size/d_sum, y_val*cell_size/d_sum
end

function covariance(d::Matrix{Float64}, L::Float64)
	mu_x, mu_y = centroid(d, L)
	c_xx = c_xy = c_yy = 0.0
	n = size(d,1)
	cell_size = L / n
	d_sum = 0.0
	for xi = 1:n
		for yi = 1:n
			x = (xi-0.5)*cell_size
			y = (yi-0.5)*cell_size

			c_xx += d[xi,yi] * x * x
			c_yy += d[xi,yi] * y * y
			c_xy += d[xi,yi] * (x - mu_x) * (y - mu_y)
			d_sum += d[xi,yi]
		end
	end
	c_xx = c_xx/d_sum - (mu_x * mu_x)
	c_yy = c_yy/d_sum - (mu_y * mu_y)
	c_xy = c_xy/d_sum

	# add 1e-3 to ensure we are positive definite
	return [c_xx+1e-3 c_xy; c_xy c_yy+1e-3]
end
