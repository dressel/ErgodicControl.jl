######################################################################
# includes stuff I don't want to need distributions for...
######################################################################
export traj2mat

function my_pdf(x::NTuple{2,Float64},m::Vector{Float64},S::Matrix{Float64})
	dx = x[1] - m[1]
	dy = x[2] - m[2]

	Si = inv(S)
	e_stuff = dx*Si[1,1]*dx + dx*Si[1,2]*dy + dy*Si[2,1]*dx + dy*Si[2,2]*dy
	return det(sqrtm(Si)) * exp(-0.5 * e_stuff) / (2.0 * pi)
end

function my_pdf(x::VF, m::VF, S::MF)
	k = length(x)
	xm = x - m
	return exp(-0.5 * dot(xm, inv(S)*xm)) / sqrt((2pi)^k * det(S))
end

# In order to be faster, precompute inverse and denominator.
# That way, I don't have to recompute them for every cell in domain.
function fast_pdf(x::VF, m::VF, inv_S::MF, den::Float64)
	xm = x - m
	return exp(-0.5 * dot(xm, inv_S*xm)) / den
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
function directional_derivative(g_f1::VVF, g_f2::VVF, u1::VVF, u2::VVF)
	N = length(g_f2)
	dd = 0.0
	for i = 1:N
		dd += dot(g_f1[i], u1[i]) + dot(g_f2[i], u2[i])
	end
	dd += dot(g_f1[N+1], u1[n+1])
	return dd
end
function directional_derivative(g_f1::Matrix{Float64}, g_f2::Matrix{Float64}, u1::VVF, u2::VVF)
	#N = length(g_f1)
	N = size(g_f2, 2)
	#n = size(u1[1])
	#m = size(u2[1])
	#println("N = ", N)
	dd = 0.0
	for i = 1:N
		# old way
		#dd += g_f1[1,i]*u1[i][1] + g_f1[2,i]*u1[i][2]
		#dd += g_f2[1,i]*u2[i][1] + g_f2[2,i]*u2[i][2]

		dd += dot(g_f1[:,i], u1[i]) + dot(g_f2[:,i], u2[i])
	end
	dd += dot(g_f1[:,N+1], u1[N+1])
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
function scaled_dd(g_f1::VVF, g_f2::VVF, u1::VVF, u2::VVF)
	dd = directional_derivative(g_f1, g_f2, u1, u2)
	norm_factor = sqrt( dot(u1,u1) + dot(u2,u2) )
	return dd / norm_factor
end

function scaled_dd(g_f1::Matrix{Float64}, g_f2::Matrix{Float64}, u1::VVF, u2::VVF)
	dd = directional_derivative(g_f1, g_f2, u1, u2)
	norm_factor = sqrt( dot(u1,u1) + dot(u2,u2) )
	return dd / norm_factor
end

function scaled_dd(g_f1::Matrix{Float64}, g_f2::Matrix{Float64}, u1::Matrix{Float64}, u2::Matrix{Float64})
	dd = directional_derivative(g_f1, g_f2, u1, u2)
	norm_factor = sqrt( sum(u1.*u1) + sum(u2.*u2) )
	return dd / norm_factor
end


# TODO: why exactly do I do this?
# "normalizes" a matrix so (bins * bins) / sum(mat) = 1.0
function normalize!(mat::Matrix{Float64}, dA::Float64)
	num_x, num_y = size(mat)
	#c = (bins * bins) / sum(mat)
	c = 1.0 / (sum(mat) * dA)	# LD 3/02/2017
	for xi = 1:num_x
		for yi = 1:num_y
			mat[xi,yi] *= c
		end
	end
end
function normalize!(mat::Array{Float64,3}, dV::Float64)
	num_x, num_y, num_z = size(mat)
	#c = (bins * bins * bins) / sum(mat)
	c = 1.0 / (sum(mat) * dV)		# LD change 3/02/2017
	for xi = 1:num_x
		for yi = 1:num_y
			for zi = 1:num_z
				mat[xi,yi,zi] *= c
			end
		end
	end
end
export normalize!

"""
`mat2traj(mat::Matrix{Float64})`

Converts a matrix to a vector of Float64 vectors (VVF's).
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

"""
`traj2mat(traj::VVF)`

Converts vector of vector of floats into `N x n` matrix, where `N` is the number of points in trajectory and `n` is the dimensionality of the state.
"""
function traj2mat(traj::VVF)
	N = length(traj)
	n = length(traj[1])
	mat = zeros(N, n)
	for i = 1:N
		mat[i,:] = traj[i]
	end
	return mat
end



export vvvf2vvf
function vvvf2vvf(xd::VVVF, ud::VVVF)
	N = length(ud[1])
	xdcat = hcat(xd...)
	udcat = hcat(ud...)
	x = VVF(N+1)
	u = VVF(N)
	for n = 1:N
		x[n] = vcat(xdcat[n,:]...)
		u[n] = vcat(udcat[n,:]...)
	end
	x[N+1] = vcat(xdcat[N+1,:]...)
	return x,u
end
