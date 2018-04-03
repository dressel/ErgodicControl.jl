# linear.jl
# testing pto_linear which should be faster
em = ErgodicManagerR2("single gaussian", K=5, bins=100)

x0 = [0.4,0.1]
N = 40
h = 0.1

tm = TrajectoryManager(x0, h, N, ConstantInitializer([0.0,0.0]))

xd, ud = pto_linear(em, tm)
