######################################################################
# smc_trajectory.jl
#
# Here I try to apply the original way of doing it.
######################################################################

export smc_trajectory

function smc_trajectory(em::ErgodicManagerR2, tm::TrajectoryManager; verbose::Bool=true, umax::Float64=1.0)


    ud = Array{Vector{Float64}}(undef, tm.N)#VVF(tm.N)
    xd = Array{Vector{Float64}}(undef, tm.N+1)#VVF(tm.N+1)
    xd[1] = deepcopy(tm.x0)
    ck = zeros(em.K+1, em.K+1)
    for n = 1:tm.N
        # compute B
        Bx, By = compute_B(em, xd, n, tm.h, ck)

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
function compute_B(em::ErgodicManagerR2, xd::VVF, n::Int, h::Float64, ck::Matrix{Float64})
    Lx = em.domain.lengths[1]
    Ly = em.domain.lengths[2]
    dxn = xd[n][1] - x_min(em)
    dyn = xd[n][2] - y_min(em)
    Bx = By = 0.0
    for k1 = 0:em.K, k2 = 0:em.K
        hk = em.hk[k1+1, k2+1]
        cx = cos(k1*pi*dxn/Lx)
        cy = cos(k2*pi*dyn/Ly)

        # compute ck
        fk = ( (n-1)*ck[k1+1,k2+1] + cx*cy/hk ) / n
        ck[k1+1, k2+1] = fk


        # t * Lambda_k * (c_k - phi_k)
        LS = em.Lambda[k1+1,k2+1] * h * n * (fk - em.phik[k1+1, k2+1])

        # multiplying by gradient
        Bx += -LS * k1*pi * sin(k1*pi*dxn/Lx) * cy / (hk*Lx)
        By += -LS * k2*pi * cx * sin(k2*pi*dyn/Ly) / (hk*Ly)
    end
    return Bx, By
end


function smc_trajectory(em::ErgodicManagerR3, tm::TrajectoryManager; verbose::Bool=true, umax::Float64=1.0)

    ud = VVF(tm.N)
    xd = VVF(tm.N+1)
    xd[1] = deepcopy(tm.x0)
    ck = zeros(em.K+1, em.K+1, em.K+3)
    for n = 1:tm.N
        # compute B
        Bx, By, Bz = compute_B(em, xd, n, tm.h, ck)

        # normalize
        den = sqrt(Bx*Bx + By*By + Bz*Bz)
        ux = -umax * Bx / den
        uy = -umax * By / den
        uz = -umax * Bz / den

        ud[n] = [ux, uy, uz]
        xd[n+1] = integrate(tm, xd[n], ud[n])
    end

    return xd, ud
end


function compute_B(em::ErgodicManagerR3, xd::VVF, n::Int, h::Float64, ck::Array{Float64,3})
    Lx = em.domain.lengths[1]
    Ly = em.domain.lengths[2]
    Lz = em.domain.lengths[3]

    dxn = xd[n][1] - x_min(em)
    dyn = xd[n][2] - y_min(em)
    dzn = xd[n][3] - z_min(em)

    Bx = By = Bz = 0.0

    for k1 = 0:em.K, k2 = 0:em.K, k3 = 0:em.K
        hk = em.hk[k1+1, k2+1, k3+1]
        cx = cos(k1*pi*dxn/Lx)
        cy = cos(k2*pi*dyn/Ly)
        cz = cos(k3*pi*dzn/Lz)

        # compute ck
        fk = ( (n-1)*ck[k1+1,k2+1,k3+1] + cx*cy*cz/hk ) / n
        ck[k1+1, k2+1, k3+1] = fk


        # t * Lambda_k * (c_k - phi_k)
        LS = em.Lambda[k1+1,k2+1,k3+1]*h*n*(fk - em.phik[k1+1, k2+1, k3+1])

        # multiplying by gradient
        Bx += -LS * k1*pi * sin(k1*pi*dxn/Lx) * cy * cz / (hk*Lx)
        By += -LS * k2*pi * cx * sin(k2*pi*dyn/Ly) * cz / (hk*Ly)
        Bz += -LS * k3*pi * cx * cy * sin(k3*pi*dzn/Lz) / (hk*Lz)
    end
    return Bx, By, Bz
end
