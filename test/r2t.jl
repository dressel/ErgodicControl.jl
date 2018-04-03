# Generate the distribution
N = 80
dt = 0.5
T = N*dt
d = Domain([1,1], [100,100])
cov = 0.010 * eye(2)
phi = zeros(100,100,N+1)
for i = 1:N+1
    mui = (.7*(i-1)/N + .15) * ones(2)
    phi[:,:,i] = gaussian(d, mui, cov)
end
ErgodicControl.normalize!(phi, d.cell_size / (N+1))

# Now let's create the ergodic manager in R3
K = 5
em = ErgodicManagerR2T(d, phi, K)

# trajectory params
x0 = [0.49, 0.01]
tm = TrajectoryManager(x0, dt, N, ConstantInitializer([0.,0.]))
tm.R = .1*eye(2)

# I call this second Armijo
tm.descender = ArmijoLineSearch(1,1e-4)

# trajectory generation and plotting
mi = 1000
ddc = 1e-5
v = true
xd,ud = pto_trajectory(em, tm, dd_crit=ddc, max_iters=mi, verbose=v)
