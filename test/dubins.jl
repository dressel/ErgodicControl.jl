em = ErgodicManagerR2("double gaussian", K=5, bins=100)

x0 = [0.5,0.01,pi/4]
N = 40
h = 0.1

tm = TrajectoryManager(x0, h, N, ConstantInitializer([0.0000]))

# things needed for dynamics
tm.dynamics = DubinsDynamics(0.3,0.1)
tm.Qn = eye(3)
tm.R = 0.01 * eye(1)
tm.Rn = 1 * eye(1)

xd, ud = pto_trajectory(em, tm)

tm.barrier_cost = 1
xd, ud = pto_trajectory(em, tm)
