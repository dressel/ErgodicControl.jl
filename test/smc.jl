em = ErgodicManagerR2("double gaussian", K=5, bins=100)

# Trajectory manager
x0 = [0.89,0.01]
N = 400
h = 0.1
ci = ConstantInitializer([-0.01,0.01])
tm = TrajectoryManager(x0, h, N, ci)

xd, ud = smc_trajectory(em, tm, umax=.15)




# now test with two agents...

# Set up different domains with different discretizations
d = Domain([1,1], 100)
num_agents = 2

# Set up distribution and ergodic manager
K = 5
means = [[.3,.7], [.7,.3]]
Sigmas = [.025*eye(2), .025*eye(2)]
phi = gaussian(d, means, Sigmas)
em = ErgodicManagerR2(d, phi, K)

# Set up first trajectory manager
x0 = [0.49,0.01]
N = 50
h = 0.6
ci = ConstantInitializer([0.0, 0.0])
tm1 = TrajectoryManager(x0, h, N, ci)
dynamics!(tm1, SingleIntegrator(2,h))

# second tm is like the first, but different starting point
tm2 = deepcopy(tm1)
tm2.x0 = [.79,.99]

# array of trajectory managers
vtm = [tm1, tm2]

# Generate the trajectories
ddc = 1e-4
xd, ud = smc_trajectory(em, vtm, umax=.15)
