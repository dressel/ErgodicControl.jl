######################################################################
# initializer.jl
# different ways to initalize a trajectory
######################################################################

abstract Initializer

"""
`ri = CornerInitializer()`

Takes a trajectory to the farthest corner.
"""
type RandomInitializer <: Initializer end


# TODO: change rand() call so we can account for various
function initialize(ri::RandomInitializer, em::ErgodicManager, x0::Vector{Float64}, h::Float64, N::Int)
	xd = [[x0[1], x0[2]] for i = 1:N+1]
	points = Array(Vector{Float64}, N)
	for i = 1:N
		points[i] = [rand(), rand()]
	end

	tsp_nn!(xd, points)

	ud = compute_controls(xd, h)

	return xd, ud
end


"""
`ci = CornerInitializer()`

Takes a trajectory to the farthest corner.
"""
type CornerInitializer <: Initializer end

# Assumes we our domain is 2d
# corner 1 = 0,0
# corner 2 = 1,0
# corner 3 = 0,1
# corner 4 = 1,1
function initialize(ci::CornerInitializer, em::ErgodicManager, x0::Vector{Float64}, h::Float64, N::Int)

	# dimensionality of state space
	n = length(x0)

	# first corner: 0,0
	dx1 = 0.0 - x0[1]
	dy1 = 0.0 - x0[2]
	dc1 = sqrt(dx1*dx1 + dy1*dy1)

	# second corner: 1,0
	dx2 = 1.0 - x0[1]
	dy2 = 0.0 - x0[2]
	dc2 = sqrt(dx2*dx2 + dy2*dy2)

	# third corner: 0,1
	dx3 = 0.0 - x0[1]
	dy3 = 1.0 - x0[2]
	dc3 = sqrt(dx3*dx3 + dy3*dy3)

	# fourth corner: 1,1
	dx4 = 1.0 - x0[1]
	dy4 = 1.0 - x0[2]
	dc4 = sqrt(dx4*dx4 + dy4*dy4)

	bi = indmax([dc1,dc2,dc3,dc4])

	bc = [0.,0.]
	if bi == 2
		bc = [1.,0.]
	elseif bi == 3
		bc = [0.,1.]
	elseif bi == 4
		bc = [1.,1.]
	end

	x_step = (bc[1] - x0[1]) / N
	y_step = (bc[2] - x0[2]) / N

	xd = Array(Vector{Float64}, N+1)
	for i = 0:N
		xd[i+1] = zeros(n)
		xd[i+1][1] = x0[1] + i*x_step
		xd[i+1][2] = x0[2] + i*y_step
	end

	ud = compute_controls(xd, h)

	return xd, ud
end

"""
`ci = ConstantInitializer()`

Just takes a constant action.
"""
type ConstantInitializer <: Initializer
	action::Vector{Float64}
end

#function initialize(ci::ConstantInitializer, em::ErgodicManager, x0::Vector{Float64}, h::Float64, N::Int)
function initialize(ci::ConstantInitializer, em::ErgodicManager, x0::Vector{Float64}, h::Float64, N::Int, A::Matrix{Float64}, B::Matrix{Float64})
	n = length(x0)
	m = length(ci.action)
	xd = Array(Vector{Float64}, N+1)
	ud = Array(Vector{Float64}, N+1)
	xd[1] = x0
	ud[1] = ci.action
	for i = 1:N
		xd[i+1] = A*xd[i] + B*ud[i]
		ud[i+1] = deepcopy(ci.action)
	end

	return xd, ud
end
