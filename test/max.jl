######################################################################
# max.jl
######################################################################
using ErgodicControl

# Ergodic Manager
d = Domain([1,1], 100)
means = [[.2,.5], [.8,.5]]
Sigmas = [.015*eye(2), .015*eye(2)]
weights = [1., 2.]
phi = gaussian(d, means, Sigmas, weights)
em = ErgodicManagerR2(d, phi, 5)

# Trajectory Manager
x0 = [0.49,0.01]
x0 = [0.49,0.01, 0, 0]
N = 30
h = 0.5
ci = ConstantInitializer([0.000, 0.000])
tm = TrajectoryManager(x0, h, N, ci)
dynamics!(tm, DoubleIntegrator(2,h))

# max trajectory
mi = 1000
tm.R = 50eye(2)
xd, ud = max_trajectory(em, tm, means, Sigmas, weights, max_iters=mi)
