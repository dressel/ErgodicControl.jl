==========================
Ergodic Control Overview
==========================
Ergodic control is concerned with the generation of ergodic trajectories.
This section defines ergodicity and shows how many of these trajectories are generated.

Ergodicity
===========
A trajectory is ergodic with respect to a distribution if its time-averaged statistics match the distribution's spatial statistics. 
In an ergodic trajectory, the time spent in a region is proportional to the distribution's density in the region.
If we have some domain :math:`X`, we denote the spatial distribution :math:`\phi`, and :math:`\phi(x)` describes the density at a point :math:`x\in X`.

A commonly used metric for ergodicity uses Fourier coefficients [1].
The spatial distribution is decomposed into Fourier coefficients :math:`\phi_k`:

.. math:: \phi_k(x) = \int_X \phi(x) F_k(x) dx,

where :math:`k` is a wave-number vector with dimensionality equal to the spatial domain's.
That is, if the spatial domain has two dimensions, then :math:`k` is a vector of length 2.

The function :math:`F_k(x)` is as follows:

.. math:: F_k(x) = \frac{1}{h_k}\prod_{i=1}^n \cos \left(\frac{k_i\pi}{L_i} x_i\right)

If we have a trajectory :math:`x(t)`, we can find the Fourier coefficients :math:`c_k` of the trajectory:

.. math:: c_k = \frac{1}{T}\int_0^T F_k(x(t))dt,

The ergodic metric :math:`\mathcal{E}` is a measure of the difference between trajectory and distribution coefficients:

.. math:: \mathcal{E} = \sum_k \Lambda_k | c_k - \phi_k |^2,

where :math:`\Lambda_k` is weighted to favor low-frequency features.

.. math:: \Lambda_k = \frac{1}{\left(1 + ||k||^2\right)^(n+1)/2},

where :math:`s = (n+1)/2` and :math:`n` is the number state variables in our distribution; for distributions over :math:`mathbb{R}^2,\ n = 2`.


Trajectory Generation
======================
We want to generate trajectories with low ergodic score, because an ergodic score of zero implies perfect ergodicity.


Bibliography
==============
[1] G. Mathew and I. Mezic, "Metrics for ergodicity and design of ergodic dynamics for multi-agent systems"
