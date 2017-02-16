######################################################################
# corner_initializer.jl
######################################################################

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
function initialize(ci::CornerInitializer, em::ErgodicManager, tm::TrajectoryManager)

	# dimensionality of state space
	x0 = tm.x0
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

	# could now do point initializer to the selected corner
	return initialize(PointInitializer(bc), em, tm)
end
