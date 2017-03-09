######################################################################
# smc_trajectory.jl
#
# Here I try to apply the original way of doing it.
######################################################################

export smc_trajectory

function smc_trajectory(em::ErgodicManagerR2, tm::TrajectoryManager; verbose::Bool=true, umax::Float64=1.0) 

	ud = VVF(tm.N)
	xd = VVF(tm.N+1)
	xd[1] = deepcopy(tm.x0)
	for n = 1:tm.N
		# compute B
		Bx, By = compute_B(em, xd, n, tm.h)

		# normalize
		den = sqrt(Bx*Bx + By*By)
		ux = -umax * Bx / den
		uy = -umax * By / den

		ud[n] = [ux, uy]
		xd[n+1] = integrate(tm, xd[n], ud[n])
	end

	return xd, ud
end

# TODO: this can be done a lot quicker
function compute_B(em::ErgodicManagerR2, xd::VVF, n::Int, h::Float64)
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

			# compute fk
			fk = 0.0
			for i = 1:n
				dx = xd[i][1] - x_min(em)
				dy = xd[i][2] - y_min(em)
				fk += cos(k1*pi*dx/Lx) * cos(k2*pi*dy/Ly)
			end
			fk /= (hk*n)
			LS = em.Lambda[k1+1,k2+1] * h * n * (fk - em.phik[k1+1, k2+1])

			Bx += -LS * k1*pi * sin(k1*pi*dxn/Lx) * cy / (hk*Lx)
			By += -LS * k2*pi * cx * sin(k2*pi*dyn/Ly) / (hk*Ly)
		end
	end
	return Bx, By
end
