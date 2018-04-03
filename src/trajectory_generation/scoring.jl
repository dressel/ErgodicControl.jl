######################################################################
# scoring.jl
# Computes ergodic, total, and control scores
######################################################################

# These are required because multi-agent is different
function ergodic_score(em::ErgodicManager, traj::VVF, d::Dynamics)
	return ergodic_score(em, traj)
end
function ergodic_score(em::ErgodicManager, traj::VVF, gd::GroupDynamics)
	xds = vvf2vvvf(traj, gd)
	ck = decompose(em, xds)
	return ergodic_score(em, ck)
end
function ergodic_score(em::ErgodicManager, traj::VVF, vtm::VTM)
	xds = vvf2vvvf(traj, vtm)
	ck = decompose(em, xds)
	return ergodic_score(em, ck)
end





# TODO: uncomment and test for time, comparing to new method in
#  ergodic_manager.jl
# It didn't matter for other types, but maybe it will here because of
#  complex numbers and stuff?
function ergodic_score(em::ErgodicManagerSE2, ck::Array{Complex{Float64},3})
	val = 0.0

    # TODO: Must be scrubbed since em.M became range and not int
	#for m = 0:em.M
	#	for n = 0:em.N
	#		for p = 0:em.P
    for (mi,m) in enumerate(em.M)
        for (ni,n) in enumerate(em.N)
            for (pl,p) in enumerate(em.P)
				d = em.phik[mi,ni,pl] - ck[mi,ni,pl]
				dr = real(d)
				di = imag(d)
				val += em.Lambda[mi,ni,pl] * (dr*dr + di*di)
			end
		end
	end
	return val
end

"""
`control_score(ud::VVF, R, h)`

Assumes only non-zero elements of `R` are corners (not anymore I think)

`control_score(tm::TrajectoryManager,u::VVF) = control_score(u, tm.R, tm.h)`
"""
function control_score(ud::VVF, R::Matrix{Float64}, h::Float64)
	cs = 0.0
	num_u = length(ud[1])
	for ui in ud

		for j = 1:num_u
		   cs += R[j,j] * ui[j] * ui[j]
		end

		# TODO: this is how I do it. Benchmark it
		#cs += dot(ui, R*ui)

		# old way I suppose. delete this dumb stuff
		#cs += R[1,1] * ui[1] * ui[1]
		#cs += R[2,2] * ui[2] * ui[2]
	end
	return 0.5 * h * cs
end

function control_score(tm::TrajectoryManager, ud::VVF)
	return control_score(ud, tm.R, tm.h)
end


# TODO: do I even use these functions?
control_score(ud::VVF) = control_score(ud, eye(2), 1.0)
function control_score(ud::VVF, N::Int)
	cs = 0.0
	for n = 1:N
		num_u = length(ud[n])
		for j = 1:num_u
			cs += ud[n][j] * ud[n][j]
		end
		#cs += ud[n][1] * ud[n][1]
		#cs += ud[n][2] * ud[n][2]
	end
	return cs
end

function control_effort(ud::VVF, N::Int=length(ud))
	cs = 0.0
	for n = 1:N
		udx = ud[n][1]
		udy = ud[n][1]
		cs += sqrt(udx*udx + udy*udy)
	end
	return cs
end


function barrier_score(em::ErgodicManager, xd::VVF, c::Float64)
	if c == 0.0; return 0.0; end

	bs = 0.0
	xmax = x_max(em)
	ymax = y_max(em)
	xmin = x_min(em)
	ymin = y_min(em)

	for xi in xd
		if (xi[1] > xmax)
			dx = xi[1] - xmax
			bs += c * dx * dx
		elseif (xi[1] < xmin)
			dx = xi[1] - xmin
			bs += c * dx * dx
		end
		if (xi[2] > ymax)
			dy = xi[2] - ymax
			bs += c * dy * dy
		elseif (xi[2] < ymin)
			dy = xi[2] - ymin
			bs += c * dy * dy
		end
	end
	return bs
end


"""
`total_score(em::ErgodicManager, tm::TrajectoryManager, xd::VVF, ud::VVF)`

Computes the total score `q*ergodic_score + sum_n h/2 un'Rn un`
"""
function total_score(em::ErgodicManager, tm::TrajectoryManager, xd::VVF, ud::VVF)
	es = tm.q * ergodic_score(em, xd, tm.dynamics)
	cs = control_score(ud, tm.R, tm.h)
	bs = barrier_score(em, xd, tm.barrier_cost)
	return es + cs + bs
end

"""
`all_scores(em::ErgodicManager, tm::TrajectoryManager, xd::VVF, ud::VVF)`
"""
function all_scores(em::ErgodicManager, tm::TrajectoryManager, xd::VVF, ud::VVF)
	es = tm.q * ergodic_score(em, xd, tm.dynamics)
	cs = control_score(ud, tm.R, tm.h)
	ts = es + cs
	return es, cs, ts
end
