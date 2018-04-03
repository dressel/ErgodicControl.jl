em = ErgodicManagerR2("double gaussian", K=5, bins=100)

# Trajectory manager
x0 = [0.89,0.01]
N = 400
h = 0.1
ci = ConstantInitializer([-0.01,0.01])
tm = TrajectoryManager(x0, h, N, ci)

xd, ud = smc_trajectory(em, tm, umax=.15)
