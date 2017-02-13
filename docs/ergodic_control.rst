==========================
Ergodic Control Overview
==========================

Brief overview of ergodic control gathering.


Ergodicity
===========
A trajectory is ergodic with respect to a distribution if its time-averaged statistics match the distribution's spatial statistics. 
In an ergodic trajectory, the time spent in a region is proportional to the distribution's density in the region.
If we have some domain :math:`X`, we denote the spatial distribution :math:`\phi`, and :math:`\phi(x)` describes the density at a point :math:`x\in X`.

A commonly used metric for ergodicity uses Fourier coefficients.
The spatial distribution is decomposed into Fourier coefficients :math:`\phi_k`:

.. math:: \phi_k(x) = \int_X \phi(x) F_k(x) dx,

where :math:`k` is a wave-number vector with dimensionality equal to the spatial domain's.
That is, if the spatial domain has two dimensions, then :math:`k` is a vector of length 2.

If we have a trajectory :math:`x(t)`, we can find the Fourier coefficients :math:`c_k` of the trajectory:

.. math:: c_k = \frac{1}{T}\int_0^T F_k(x(t))dt,

where `\Lambda_k` is weighted to favor low-frequency features.

The ergodic metric :math:`\mathcal{E}` is a measure of the difference between trajectory and distribution coefficients:

.. math:: \mathcal{E} = \sum_K \Lambda_k \abs{c_k - \phi_k}^2.


Trajectory Generation
======================
