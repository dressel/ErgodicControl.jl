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


x_min(em::ErgodicManager) = em.domain.mins[1]
y_min(em::ErgodicManager) = em.domain.mins[2]
z_min(em::ErgodicManager) = em.domain.mins[3]

x_max(em::ErgodicManager) = em.domain.maxes[1]
y_max(em::ErgodicManager) = em.domain.maxes[2]
z_max(em::ErgodicManager) = em.domain.maxes[3]

x_size(em::ErgodicManager) = em.domain.cell_lengths[1]
y_size(em::ErgodicManager) = em.domain.cell_lengths[2]
z_size(em::ErgodicManager) = em.domain.cell_lengths[3]

x_cells(em::ErgodicManager) = em.domain.cells[1]
y_cells(em::ErgodicManager) = em.domain.cells[2]
z_cells(em::ErgodicManager) = em.domain.cells[3]
