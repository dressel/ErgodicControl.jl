# domain, distribution, and ergodic manager
d = Domain([1,1,1], 100)
means = [[.2,.2,.2], [.8,.8,.2], [.5,.5,.8]]
covs = [0.01*eye(3), 0.01*eye(3), .01*eye(3)]
phi = gaussian(d, means, covs)
K = 5
em = ErgodicManagerR3(d, phi, K)

# trajectory params
x0 = [0.49, 0.01, 0.01]
dt = 0.5
N = 80
tm = TrajectoryManager(x0, dt, N, ConstantInitializer([0.0,0.0,0.0]))
dynamics!(tm, SingleIntegrator(3,dt))
tm.descender = ArmijoLineSearch(1,1e-4)

# trajectory generation and plotting
xd,ud = pto_trajectory(em, tm, dd_crit=1e-4, max_iters=1000)
