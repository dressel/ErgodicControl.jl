######################################################################
# ergodic_manager.jl
######################################################################
abstract ErgodicManager

decompose!(em::ErgodicManager) = decompose!(em, em.phi)

"""
`phi = reconstruct(em)`

Reconstructs distribution from Fourier coefficients `em.phik`.

`phi = reconstruct(em, ck::Matrix{Float64})`

Reconstructs distribution from Fourier coefficients `ck`.

Both options return a distribution `phi`, a `em.bins` by `em.bins` matrix.
The value at `phi[xi,yi]` gives the value at the cell index `xi,yi`.

`phi = reconstruct(em, xd::VV_F, start_idx=1)`

Reconstructs from trajectory `xd`. If `start_idx` is 1, then the right Riemann sum is used. If it is 0, the left Riemann sum is used.
"""
reconstruct(em::ErgodicManager) = reconstruct(em, em.phik)
