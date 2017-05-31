######################################################################
# smc_multi.jl
#
# Like smc.jl, but for multi-agent trajectories
# 
# In theory, this shouldn't be so different from single version
######################################################################

function smc_trajectory(em::ErgodicManagerR2, vtm::VTM; verbose::Bool=true, umax::Float64=1.0) 

	num_agents = length(vtm)
	N = vtm[1].N
	dt = vtm[1].h
	#ud = VVF(tm.N)
	#xd = VVF(tm.N+1)
	xds = VVVF(num_agents)
	uds = VVVF(num_agents)
	for j = 1:num_agents
		xds[j] = VVF(N+1)
		uds[j] = VVF(N)
		xds[j][1] = deepcopy(vtm[j].x0)
	end

	ck = zeros(em.K+1, em.K+1)
	for n = 1:N

		# update ck
		update_ck(em, xds, n, ck)

		# compute B
		for j = 1:num_agents
			Bx, By = compute_Bj(em, xds[j], n, dt, ck, num_agents)

			# normalize
			den = sqrt(Bx*Bx + By*By)
			ux = -umax * Bx / den
			uy = -umax * By / den

			uds[j][n] = [ux, uy]
			xds[j][n+1] = integrate(vtm[j], xds[j][n], uds[j][n])
		end
	end

	return vvvf2vvf(xds, uds)
end

function update_ck(em::ErgodicManagerR2, xds::VVVF, n::Int, ck::MF)
	num_agents = length(xds)
	Lx = em.domain.lengths[1]
	Ly = em.domain.lengths[2]
	for k1 = 0:em.K
		for k2 = 0:em.K
			hk = em.hk[k1+1, k2+1]
			ck[k1+1, k2+1] *= (n-1)
			for j = 1:num_agents
				dxn = xds[j][n][1] - x_min(em)
				dyn = xds[j][n][2] - y_min(em)

				cx = cos(k1*pi*dxn/Lx)
				cy = cos(k2*pi*dyn/Ly)

				# compute ck
				ck[k1+1, k2+1] += cx*cy/(hk*num_agents)
			end
			ck[k1+1, k2+1] /= n
		end
	end
end

# TODO: this can be done a lot quicker
function compute_Bj(em::ErgodicManagerR2, xd::VVF,n::Int,h::Float64,ck::MF, num_agents::Int)
	Lx = em.domain.lengths[1]
	Ly = em.domain.lengths[2]
	dxn = xd[n][1] - x_min(em)
	dyn = xd[n][2] - y_min(em)
	Bx = By = 0.0
	for k1 = 0:em.K
		for k2 = 0:em.K
			hk = em.hk[k1+1, k2+1]
			cx = cos(k1*pi*dxn/Lx)
			cy = cos(k2*pi*dyn/Ly)

			# compute ck
			#fk = ( (n-1)*ck[k1+1,k2+1] + cx*cy/hk ) / n
			#ck[k1+1, k2+1] = fk


			# t * Lambda_k * (c_k - phi_k)
			de = ck[k1+1, k2+1] - em.phik[k1+1,k2+1]
			LS = em.Lambda[k1+1,k2+1] * h*n * num_agents * de

			# multiplying by gradient
			Bx += -LS * k1*pi * sin(k1*pi*dxn/Lx) * cy / (hk*Lx)
			By += -LS * k2*pi * cx * sin(k2*pi*dyn/Ly) / (hk*Ly)
		end
	end
	return Bx, By
end
