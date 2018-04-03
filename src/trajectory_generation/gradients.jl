######################################################################
# gradients.jl
# generation of gradients
#
# TODO: this whole file needs work
#  While technically correct, it is an eyesore and can probably
#   be greatly simplified.
######################################################################

# Only first two states matter for ergodic score and barrier penalty
# Assumes ad has been initialized with zeros; that is, ad[3:end, ni] = 0.0
function gradients!(ad::MF, bd::MF, em::ErgodicManager, tm::TrajectoryManager, xd::VVF, ud::VVF)
	ck = decompose(em, xd)
	gradients!(ad, bd, em, tm, xd, ud, ck)
end


function gradients!(ad::MF, bd::MF, em::ErgodicManager, tm::TrajectoryManager, xd::VVF, ud::VVF, ck)
	ni =  1
	#ck = decompose(em, xd)
	for n = 0:(tm.N-1)

		# ergodic gradients
		an = compute_ans(em, xd, tm, n, ck)
		for i = 1:length(an)
			ad[i,ni] = an[i]
		end

		# quadratic boundary
		if tm.barrier_cost > 0.0
			xnx = xd[n+1][1]
			xny = xd[n+1][2]
			xmax = x_max(em)
			xmin = x_min(em)
			ymax = y_max(em)
			ymin = y_min(em)
			if (xnx > xmax)
				ad[1,ni] += tm.barrier_cost * 2.0 * (xnx - xmax)
			elseif xnx < xmin
				ad[1,ni] += tm.barrier_cost * 2.0 * (xnx - xmin)
			end
			if xny > ymax
				ad[2,ni] += tm.barrier_cost * 2.0 * (xny - ymax)
			elseif xny < ymin
				ad[2,ni] += tm.barrier_cost * 2.0 * (xny - ymin)
			end
		end

		# control gradients
		bd[:,ni] = tm.h * tm.R * ud[ni]

		ni += 1
	end

	n = tm.N
	an = compute_ans(em, xd, tm, n, ck)
	for i = 1:length(an)
		ad[i,ni] = an[i]
	end
	if tm.barrier_cost > 0.0
		xnx = xd[n+1][1]
		xny = xd[n+1][2]
		xmax = x_max(em)
		xmin = x_min(em)
		ymax = y_max(em)
		ymin = y_min(em)
		if (xnx > xmax)
			ad[1,ni] += tm.barrier_cost * 2.0 * (xnx - xmax)
		elseif xnx < xmin
			ad[1,ni] += tm.barrier_cost * 2.0 * (xnx - xmin)
		end
		if xny > ymax
			ad[2,ni] += tm.barrier_cost * 2.0 * (xny - ymax)
		elseif xny < ymin
			ad[2,ni] += tm.barrier_cost * 2.0 * (xny - ymin)
		end
	end
end


# This version is for multi-agent method
function gradients!(ad::MF, bd::MF, em::ErgodicManager, vtm::Vector{TrajectoryManager}, xd::VVF, ud::VVF)

	# must first create xds, where each agent trajectory is its own vector
	xds, uds = vvf2vvvf(xd, ud, vtm)

	ck = decompose(em, xds)
	n_idx = 1
	m_idx = 1
	N = length(xd) - 1
	num_agents = length(vtm)
	for j = 1:num_agents
		tm = vtm[j]
		adt = zeros(tm.dynamics.n, N+1)
		bdt = zeros(tm.dynamics.m, N)
		n_rows = tm.dynamics.n
		m_rows = tm.dynamics.m
		gradients!(adt, bdt, em, tm, xds[j], uds[j], ck)
		ad[n_idx:(n_idx+n_rows-1), :] = adt / num_agents
		bd[m_idx:(m_idx+m_rows-1), :] = bdt / num_agents
		n_idx += n_rows
		m_idx += m_rows
	end
end


# TODO: make this functional
function add_barrier_gradient!(ad::MF, em::ErgodicManager, tm::TrajectoryManager)
	if tm.barrier_cost > 0.0
		xnx = xd[n+1][1]
		xny = xd[n+1][2]
		xmax = x_max(em)
		xmin = x_min(em)
		ymax = y_max(em)
		ymin = y_min(em)
		if (xnx > xmax)
			ad[1,ni] += tm.barrier_cost * 2.0 * (xnx - xmax)
		elseif xnx < xmin
			ad[1,ni] += tm.barrier_cost * 2.0 * (xnx - xmin)
		end
		if xny > ymax
			ad[2,ni] += tm.barrier_cost * 2.0 * (xny - ymax)
		elseif xny < ymin
			ad[2,ni] += tm.barrier_cost * 2.0 * (xny - ymin)
		end
	end
end

# computes a_n, derivative of c_k wrt x_n
# returns tuple containing elements of an
# TODO: I shouldn't need to use tm
function compute_ans(em::ErgodicManagerR2, xd::VVF, tm::TrajectoryManager, n::Int, ck::Matrix{Float64})
	x = xd[n + 1][1]
	y = xd[n + 1][2]

	Lx = em.domain.lengths[1]
	Ly = em.domain.lengths[2]

	an_x = 0.0
	an_y = 0.0
	 
	xm = x_min(em)
	ym = y_min(em)

	for k1 = 0:em.K
		for k2 = 0:em.K
			hk = em.hk[k1+1,k2+1]

			dFk_dxn1 = -k1*pi*sin(k1*pi*(x-xm)/Lx)*cos(k2*pi*(y-ym)/Ly) / (hk*Lx)
			dFk_dxn2 = -k2*pi*cos(k1*pi*(x-xm)/Lx)*sin(k2*pi*(y-ym)/Ly) / (hk*Ly)

			c = em.Lambda[k1+1,k2+1] * (ck[k1+1,k2+1] - em.phik[k1+1,k2+1])
			an_x += c*dFk_dxn1
			an_y += c*dFk_dxn2
		end
	end
	an_x *= 2.0/(tm.N+1)
	an_y *= 2.0/(tm.N+1)
	return an_x, an_y
end

function compute_ans(em::ErgodicManagerR3, xd::VVF, tm::TrajectoryManager, n::Int, ck::Array{Float64,3})
	x = xd[n + 1][1]
	y = xd[n + 1][2]
	z = xd[n + 1][3]

	Lx = em.domain.lengths[1]
	Ly = em.domain.lengths[2]
	Lz = em.domain.lengths[3]

	an_x = 0.0
	an_y = 0.0
	an_z = 0.0
	 
	xm = x_min(em)
	ym = y_min(em)
	zm = z_min(em)

	for k1 = 0:em.K
		for k2 = 0:em.K
			for k3 = 0:em.K
				hk = em.hk[k1+1,k2+1,k3+1]

				c1 = cos(k1*pi*(x-xm)/Lx)
				c2 = cos(k2*pi*(y-ym)/Ly)
				c3 = cos(k3*pi*(z-zm)/Lz)

				dFk_dxn1 = -k1*pi* sin(k1*pi*(x-xm)/Lx)*c2*c3 / (hk*Lx)
				dFk_dxn2 = -k2*pi*c1* sin(k2*pi*(y-ym)/Ly)*c3 / (hk*Ly)
				dFk_dxn3 = -k3*pi*c1*c2* sin(k3*pi*(z-zm)/Lz) / (hk*Lz)

				c = em.Lambda[k1+1,k2+1,k3+1] * (ck[k1+1,k2+1,k3+1] - em.phik[k1+1,k2+1,k3+1])
				an_x += c*dFk_dxn1
				an_y += c*dFk_dxn2
				an_z += c*dFk_dxn3
			end
		end
	end
	an_x *= 2.0/(tm.N+1)
	an_y *= 2.0/(tm.N+1)
	an_z *= 2.0/(tm.N+1)
	return an_x, an_y, an_z
end


function compute_ans(em::ErgodicManagerR2T, xd::VVF, tm::TrajectoryManager, n::Int, ck::Array{Float64,3})
	x = xd[n + 1][1]
	y = xd[n + 1][2]
	t = n*tm.h

	Lx = em.domain.lengths[1]
	Ly = em.domain.lengths[2]
	T = tm.N*tm.h

	an_x = 0.0
	an_y = 0.0
	 
	xm = x_min(em)
	ym = y_min(em)

	for k1 = 0:em.K
		for k2 = 0:em.K
			for k3 = 0:em.K
				hk = em.hk[k1+1,k2+1,k3+1]

				c1 = cos(k1*pi*(x-xm)/Ly)
				c2 = cos(k2*pi*(y-ym)/Ly)
				c3 = cos(k3*pi*t/T)

				dFk_dxn1 = -k1*pi* sin(k1*pi*(x-xm)/Lx)*c2*c3 / (hk*Lx)
				dFk_dxn2 = -k2*pi*c1* sin(k2*pi*(y-ym)/Ly)*c3 / (hk*Ly)

				c = em.Lambda[k1+1,k2+1,k3+1] * (ck[k1+1,k2+1,k3+1] - em.phik[k1+1,k2+1,k3+1])
				an_x += c*dFk_dxn1
				an_y += c*dFk_dxn2
			end
		end
	end
	an_x *= 2.0/(tm.N+1)
	an_y *= 2.0/(tm.N+1)
	return an_x, an_y
end

function compute_ans(em::ErgodicManagerSE2, xd::VVF, tm::TrajectoryManager, n::Int, ck)
	x = xd[n + 1][1]
	y = xd[n + 1][2]
	z = xd[n + 1][3]

	an_x = 0.0
	an_y = 0.0
	an_z = 0.0

	i = float(im)

	# polar coords
	r2 = x*x + y*y	# range squared
	r = sqrt(r2)
	psi = atan2(y, x)
	 
    # TODO: this needs to be scrubbed real well, especially since I changed
    #        from em.M being an integer to a range
	#for m = 0:em.M
    #    for n = 0:em.N
    for (mi,m) in enumerate(em.M)
        for (ni,n) in enumerate(em.N)
			# commonly used values
			inm = i^(n-m)
			expi = exp(i * (m*psi + (n-m)*z) )
			#for p = 0:em.P
            for (pl,p) in enumerate(em.P)

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


				#c = em.Lambda[m+1,n+1,p+1] * (ck[m+1,n+1,p+1] - em.phik[m+1,n+1,p+1])
				c = em.Lambda[mi,ni,pl] * (ck[mi,ni,pl] - em.phik[mi,ni,pl])
				c_r = real(c)
				c_i = imag(c)

				# recall that c is a complex number...
				an_x += c_r*real(dF_dx) + c_i*imag(dF_dx)
				an_y += c_r*real(dF_dy) + c_i*imag(dF_dy)
				an_z += c_r*real(dF_dz) + c_i*imag(dF_dz)
			end
		end
	end
	an_x *= 2.0 / (tm.N+1)
	an_y *= 2.0 / (tm.N+1)
	an_z *= 2.0 / (tm.N+1)
	return an_x, an_y, an_z
end
