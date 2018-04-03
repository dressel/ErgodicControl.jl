# not generating a trajectory yet, let's just test ergodic manager

n = 28
K = 17

d = Domain([0,0,0], [1,1,2*pi], [n,n,36])
em = ErgodicManagerSE2(d, K)

phi = gaussian(d, [0.5,0.5,pi], 0.03*eye(3))
copy!(em.phi, phi)
normalize!(em)
decompose!(em)
