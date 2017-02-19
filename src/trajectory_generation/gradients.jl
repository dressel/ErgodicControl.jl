######################################################################
# gradients.jl
# generation of gradients
######################################################################

# Only first two states matter for ergodic score and barrier penalty
# Assumes ad has been initialized with zeros; that is, ad[3:end, ni] = 0.0
function gradients!(ad::Matrix{Float64}, bd::Matrix{Float64}, em::ErgodicManager, tm::TrajectoryManager, xd::VV_F, ud::VV_F, start_idx::Int)
	ni =  1
	for n = 0:(tm.N-1)

		# ergodic gradients
		an_x, an_y = compute_ans(em, xd, tm, n, start_idx)
		ad[1,ni] = an_x
		ad[2,ni] = an_y

		# quadratic boundary
		xnx = xd[n+start_idx+1][1]
		xny = xd[n+start_idx+1][2]
		c = 10.
		if xnx > 1
			ad[1,ni] += c * (2.0*xnx - 1.0)
		elseif xnx < 0
			ad[1,ni] += c * 2.0 * xnx
		end
		if xny > 1
			ad[2,ni] += c * (2.0*xny - 1.0)
		elseif xny < 0
			ad[2,ni] += c * 2.0 * xny
		end
		#ad[1,ni] += c * (2.0*xnx - 1.0)
		#ad[2,ni] += c * (2.0*xny - 1.0)

		# control gradients
		bd[:,ni] = tm.h * tm.R * ud[ni]
		# The code below reduces some allocation but is barbaric
		# The difference is minimal so I just use the line above
		#if tm.dynamics.m == 2
		#	bd[1,ni] = tm.h * (tm.R[1,1]*ud[ni][1] + tm.R[1,2]*ud[ni][2])
		#	bd[2,ni] = tm.h * (tm.R[2,1]*ud[ni][1] + tm.R[2,2]*ud[ni][2])
		#end
		#if tm.dynamics.m == 1
		#	bd[1,ni] = tm.h * (tm.R[1,1]*ud[ni][1])
		#end

		ni += 1
	end
end

function compute_ans(em::ErgodicManager, xd::VV_F, tm::TrajectoryManager, n::Int, start_idx::Int)
	xnx = xd[n + start_idx + 1][1]
	xny = xd[n + start_idx + 1][2]
	L = em.L

	an_x = 0.0
	an_y = 0.0
	 
	for k1 = 0:em.K
		for k2 = 0:em.K
			hk = em.hk[k1+1,k2+1]

			dFk_dxn1 = -k1*pi*sin(k1*pi*xnx/L)*cos(k2*pi*xny/L) / (hk*L)
			dFk_dxn2 = -k2*pi*cos(k1*pi*xnx/L)*sin(k2*pi*xny/L) / (hk*L)

			fk = 0.0
			for i in 0:(tm.N-1)
				x = xd[i + start_idx + 1]
				fk += cos(k1*pi*x[1]/L) * cos(k2*pi*x[2]/L) / hk
			end
			c = em.Lambdak[k1+1,k2+1] * (tm.h*fk/tm.T - em.phik[k1+1,k2+1])
			an_x += c*dFk_dxn1
			an_y += c*dFk_dxn2
		end
	end
	an_x *= 2.0*tm.h
	an_y *= 2.0*tm.h
	return an_x, an_y
end
