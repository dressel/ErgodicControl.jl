######################################################################
# gradients.jl
# generation of gradients
######################################################################

# Only first two states matter for ergodic score and barrier penalty
# Assumes ad has been initialized with zeros; that is, ad[3:end, ni] = 0.0
function gradients!(ad::Matrix{Float64}, bd::Matrix{Float64}, em::ErgodicManager, tm::TrajectoryManager, xd::VVF, ud::VVF, start_idx::Int)
	ni =  1
	ck = decompose(em, xd, start_idx)
	for n = 0:(tm.N-1)

		# ergodic gradients
		an = compute_ans(em, xd, tm, n, start_idx, ck)
		for i = 1:length(an)
			ad[i,ni] = an[i]
		end
		#an_x, an_y = compute_ans(em, xd, tm, n, start_idx)
		#ad[1,ni] = an_x
		#ad[2,ni] = an_y

		# quadratic boundary
		#if tm.barrier_cost > 0.0
		#	xnx = xd[n+start_idx+1][1]
		#	xny = xd[n+start_idx+1][2]
		#	if (xnx > em.L) || (xnx < 0.0)
		#		ad[1,ni] += tm.barrier_cost * (2.0*xnx - em.L)
		#	end
		#	if (xny > em.L) || (xny < 0.0)
		#		ad[2,ni] += tm.barrier_cost * (2.0*xny - em.L)
		#	end
		#end
		# alternate way...
		if tm.barrier_cost > 0.0
			xnx = xd[n+start_idx+1][1]
			xny = xd[n+start_idx+1][2]
			if (xnx > em.L)
				ad[1,ni] += tm.barrier_cost * (2.0*xnx - 2.0*em.L)
			elseif xnx < 0.0
				ad[1,ni] += tm.barrier_cost * 2.0*xnx
			end
			if xny > em.L
				ad[2,ni] += tm.barrier_cost * (2.0*xny - 2.0*em.L)
			elseif xny < 0.0
				ad[2,ni] += tm.barrier_cost * (2.0*xny)
			end
		end

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

function compute_ans(em::ErgodicManagerR2, xd::VV_F, tm::TrajectoryManager, n::Int, start_idx::Int, ck::Matrix{Float64})
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

			c = em.Lambda[k1+1,k2+1] * (ck[k1+1,k2+1] - em.phik[k1+1,k2+1])
			an_x += c*dFk_dxn1
			an_y += c*dFk_dxn2
		end
	end
	an_x *= 2.0*tm.h
	an_y *= 2.0*tm.h
	return an_x, an_y
end

function compute_ans(em::ErgodicManagerSE2, xd::VV_F, tm::TrajectoryManager, n::Int, start_idx::Int, ck)
	x = xd[n + start_idx + 1][1]
	y = xd[n + start_idx + 1][2]
	z = xd[n + start_idx + 1][3]

	an_x = 0.0
	an_y = 0.0
	an_z = 0.0

	i = float(im)

	# polar coords
	r2 = x*x + y*y	# range squared
	r = sqrt(r2)
	psi = atan2(y, x)
	 
	for m = 0:em.M
		for n = 0:em.N
			# commonly used values
			inm = i^(n-m)
			expi = exp(i * (m*psi + (n-m)*z) )
			for p = 0:em.P

				pr = p*r	# commonly used
				bjmn = besselj(m-n, pr)
				dJ_dpr = 0.5 * (besselj(m-n-1, pr) - besselj(m-n+1, pr))

				# TODO: we can reuse even more of this stuff...
				dF_dx = expi * dJ_dpr * p * (x/r)
				dF_dx += bjmn * expi*m*i * (-y/r2)
				dF_dx = inm * dF_dx

				dF_dy = expi * dJ_dpr * p * (y/r)
				dF_dy += bjmn * expi*m*i * (x/r2)
				dF_dy = inm * dF_dy

				dF_dz = inm * bjmn * expi * i * (n-m)


				c = em.Lambda[m+1,n+1,p+1] * (ck[m+1,n+1,p+1] - em.phik[m+1,n+1,p+1])
				c_r = real(c)
				c_i = imag(c)

				# recall that c is a complex number...
				an_x += c_r*real(dF_dx) + c_i*imag(dF_dx)
				an_y += c_r*real(dF_dy) + c_i*imag(dF_dy)
				an_z += c_r*real(dF_dz) + c_i*imag(dF_dz)
			end
		end
	end
	an_x *= 2.0 * (tm.h / tm.T)
	an_y *= 2.0 * (tm.h / tm.T)
	an_z *= 2.0 * (tm.h / tm.T)
	return an_x, an_y, an_z
end
