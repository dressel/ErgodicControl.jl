# not generating a trajectory yet, let's just test ergodic manager

n = 28
K = 17

d = Domain([0,0,0], [1,1,2*pi], [n,n,36])
em = ErgodicManagerSE2(d, K)

phi = gaussian(d, [0.5,0.5,pi], 0.03*eye(3))
copy!(em.phi, phi)
normalize!(em)
decompose!(em)

phi2 = reconstruct(em, em.phik)

x0 = [0.49, 0.01, 0.01]
dt = 0.5
N = 80
tm = TrajectoryManager(x0, dt, N, ConstantInitializer([0.0,0.0,0.0]))
dynamics!(tm, SingleIntegrator(3,dt))

xd, ud = pto_trajectory(em, tm, max_iters=20)
