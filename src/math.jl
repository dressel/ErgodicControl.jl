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

# unscaled
# g_f1 = grad_1 f
# g_f2 = grad_2 f
# u1, u2 are components of direction
# basically,
# [g_f1, g_f2]' * [u1, u2]'
function directional_derivative(g_f1::VV_F, g_f2::VV_F, u1::VV_F, u2::VV_F)
	N = length(g_f1)
	dd = 0.0
	for i = 1:N
		dd += dot(g_f1[i], u1[i]) + dot(g_f2[i], u2[i])
	end
	return dd
end
function directional_derivative(g_f1::Matrix{Float64}, g_f2::Matrix{Float64}, u1::VV_F, u2::VV_F)
	#N = length(g_f1)
	N = size(g_f1, 2)
	#println("N = ", N)
	dd = 0.0
	for i = 1:N
		#dd += dot(g_f1[i], u1[i]) + dot(g_f2[i], u2[i])
		dd += g_f1[1,i]*u1[i][1] + g_f1[2,i]*u1[i][2]
		dd += g_f2[1,i]*u2[i][1] + g_f2[2,i]*u2[i][2]
	end
	return dd
end
# TODO: turn this one into somepin correct for matrix{Float64} version
function directional_derivative(g_f1::Matrix{Float64}, g_f2::Matrix{Float64}, u1::Matrix{Float64}, u2::Matrix{Float64})
	#N = length(g_f1)
	N = size(g_f1, 2)
	dd = 0.0
	for i = 1:N
		#dd += dot(g_f1[i], u1[i]) + dot(g_f2[i], u2[i])
		dd += g_f1[1,i]*u1[1,i] + g_f1[2,i]*u1[2,i]
		dd += g_f2[1,i]*u2[1,i] + g_f2[2,i]*u2[2,i]
	end
	return dd
end
# like above, but scales the direction first
function scaled_dd(g_f1::VV_F, g_f2::VV_F, u1::VV_F, u2::VV_F)
	dd = directional_derivative(g_f1, g_f2, u1, u2)
	norm_factor = sqrt( dot(u1,u1) + dot(u2,u2) )
	return dd / norm_factor
end

function scaled_dd(g_f1::Matrix{Float64}, g_f2::Matrix{Float64}, u1::VV_F, u2::VV_F)
	dd = directional_derivative(g_f1, g_f2, u1, u2)
	norm_factor = sqrt( dot(u1,u1) + dot(u2,u2) )
	return dd / norm_factor
end

function scaled_dd(g_f1::Matrix{Float64}, g_f2::Matrix{Float64}, u1::Matrix{Float64}, u2::Matrix{Float64})
	dd = directional_derivative(g_f1, g_f2, u1, u2)
	norm_factor = sqrt( sum(u1.*u1) + sum(u2.*u2) )
	return dd / norm_factor
end


# normalizes zd and vd
# recall that zd and vd are together a direction, we normalize em both
function normalizer!(zd::VV_F, vd::VV_F)
	norm_factor = sqrt(dot(zd,zd) + dot(vd,vd))
	for i = 1:length(zd)
		zd[i][1] /= norm_factor
		zd[i][2] /= norm_factor
		vd[i][1] /= norm_factor
		vd[i][2] /= norm_factor
	end
end

# "normalizes" a matrix so (bins * bins) / sum(mat) = 1.0
function normalize!(mat::Matrix{Float64})
	bins, rar = size(mat)
	c = (bins * bins) / sum(mat)
	for i = 1:bins
		for j = 1:bins
			mat[i,j] *= c
		end
	end
end
export normalize!

"""
`mat2traj(mat::Matrix{Float64})`

Converts a matrix to a vector of Float64 vectors (VV_F's).
Assumes the matrix is `N x n`, where `n` is the dimensionality of the state space and `N` is the number of points in the trajectory.
"""
function mat2traj(mat::Matrix{Float64})
	N,n = size(mat)
	traj = Array(Vector{Float64}, N)
	for i = 1:N
		traj[i] = vec(mat[i,:])
	end
	return traj
end
