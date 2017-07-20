=========================
Ergodic Manager
=========================

The abstract :code:`ErgodicManager` type contains information about the spatial distribution. Currently, there are two sub-types. The type :code:`ErgodicManagerR2` manages distributions over :math:`\mathbb{R}^2`, and the type :code:`ErgodicManagerSE2` manages distributions over the special Euclidean group SE(2).

Fields
=========
The :code:`ErgodicManagerR2` type has the following fields:
::

	domain::Domain              # spatial domain
	K::Int                      # number of Fourier coefficients
	phi::Matrix{Float64}        # spatial distribution
	phik::Matrix{Float64}       # distribution's Fourier coefficients

	# constant regardless of phi (depend on k1,k2)
	Lambda::Matrix{Float64}
	hk::Matrix{Float64}

	# to speed up computation
	kpixl::Matrix{Float64}
	kpiyl::Matrix{Float64}


Construction
=============
An ergodic manager for :math:`\mathbb{R}^2` can be constructed with a :code:`Domain`, a distribution over that domain, and the number of Fourier coefficients.
::

    ErgodicManagerR2(d::Domain, phi::Matrix{Float64}, K::Int)


Updating Spatial Distribution
==============================
The :code:`decompose!` function decomposes an ergodic manager's spatial distribution :code:`phi` into Fourier coefficients, updating the managers :code:`phik` field:
::

    decompose!(em::ErgodicManager)


Reconstructing Spatial Distributions
=====================================
Sometimes we want to reconstruct a spatial distribution from the Fourier coefficients, to see how well the Fourier coefficients capture the distribution.
::

    phi = reconstruct(em::ErgodicManager)

If you have your own set of coefficients :code:`ck`, you can use that instead of :code:`em.phik`:
::

    phi = reconstruct(em::ErgodicManager, ck::Matrix{Float64})

You could also pass in a trajectory and :code:`reconstruct` will take care of decomposing it into coefficients and using these to generate a spatial distribution.
::

    phi = reconstruct(em::ErgodicManager, xd::VVF)



Example Managers
=================
I provide a number of pre-made ergodic managers that correspond to frequently used example domains/distributions.
::

    ErgodicManagerR2(example_name::String; K::Int=5, bins::Int=100)

Currently, there are two valid values for example_name: "single gaussian" and "double gaussian". For exmaple, you could run:
::

    em = ErgodicManagerR2("single gaussian", K=5, bin=100)

