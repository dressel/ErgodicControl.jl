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

.. math:: F_k(x) = \frac{1}{h_k}\prod_{i=1}^n \cos \left(\frac{k_i\pi}{L_i} x_i\right),

where :math:`h_k` is a normalizing factor and :math:`L_i` is the length of dimension :math:`i`.

.. math:: h_k = \left(\int_0^{L_1} \int_0^{L_2} \cos^2(\frac{k_1\pi x_1}{L_1}) \cos^2(\frac{k_2\pi x_2}{L_2})dx_1 dx_2 \right)^{1/2}

If we have a trajectory :math:`x(t)`, we can find the Fourier coefficients :math:`c_k` of the trajectory:

.. math:: c_k = \frac{1}{T}\int_0^T F_k(x(t))dt,

The ergodic metric :math:`\mathcal{E}` is a measure of the difference between trajectory and distribution coefficients:

.. math:: \mathcal{E} = \sum_k \Lambda_k | c_k - \phi_k |^2,

where :math:`\Lambda_k` is weighted to favor low-frequency features. It takes the form

.. math:: \Lambda_k = \frac{1}{\left(1 + ||k||^2\right)^{(n+1)/2}},

where :math:`n` is the number state variables in our distribution; for distributions over :math:`\mathbb{R}^2,\ n = 2`.


Ergodicity in SE(2)
=====================
A robot in SE(2) has a state :math:`x\in\mathbb{R}^3`, where :math:`x` consists of the robot's :math:`x`-location, :math:`y`-location, and heading :math:`\theta`. We might denote this state :math:`x(t) = [x_r, y_r, \theta]`.

The basis functions are as follows:

.. math:: F_{m,n,p}(x) = i^{n-m}\exp\left( i\left[m\psi + (n-m)\theta\right]\right) J_{m-n}(pr),

where :math:`(m,n,p)` are the indices along each direction, :math:`J_{m-n}` is the :math:`m-n\text{th}` order Bessel function and :math:`(r, \psi, \theta)` are the robot's state expressed in polar coordinates:

.. math:: r = \sqrt{x_r^2 + y_r^2},
.. math:: \psi = \arctan(y_r / x_r),
.. math:: \theta = \theta.

The ergodic metric is then expressed

.. math:: \mathcal{E} = \sum_{m,n,p=0}^{M,N,P} \Lambda_{m,n,p} || c_{m,n,p} - \phi_{m,n,p} ||^2,

where :math:`\Lambda_{m,n,p}` is equivalent to :math:`\Lambda_k` where the vector :math:`k=[m,n,p]`.


Trajectory Generation
======================
We want to generate trajectories with low ergodic score, because an ergodic score of zero implies perfect ergodicity.


Bibliography
==============
[1] G. Mathew and I. Mezic, "Metrics for ergodicity and design of ergodic dynamics for multi-agent systems"
