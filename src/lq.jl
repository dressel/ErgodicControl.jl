######################################################################
# lq.jl
#
# Solves LQ descent direction problems
######################################################################
export LQ, apply_LQ_gains

function LQ(A::MF, B::MF, a::MF, b::MF, Q::MF, R::MF, N::Int)
	# Also needed for LQR
	P = Array(Matrix{Float64}, N+1)
	G = Array(Matrix{Float64}, N)
	K = Array(Matrix{Float64}, N)
	P[N+1] = Q

	# Solely for LQ
	r = Array(Vector{Float64}, N+1)
	r[N+1] = zeros(2)
	C = Array(Vector{Float64}, N)

	# Sweep from the back
	for n = (N-1):-1:0
		G[n+1] = R + (B' * P[n+1+1] * B)
		K[n+1] = inv(G[n+1]) * B * P[n+1+1] * A
		P[n+1] = Q + (A' * P[n+1+1] * A) - (K[n+1]' * G[n+1] * K[n+1])
		r[n+1] = (A'-K[n+1]'*B')r[n+1+1] + .5*a[:,n+1] - .5*K[n+1]*b[:,n+1]
		C[n+1] = inv(G[n+1]) * (B'*r[n+1+1] + 0.5*b[:,n+1])
	end

	return K, C
end

function LQ(tm::TrajectoryManager, ad::MF, bd::MF) 
	LQ(tm.dynamics.A, tm.dynamics.B, ad, bd, tm.Q, tm.R, tm.N)
end

function apply_LQ_gains(A::MF, B::MF, K::Vector{MF}, C::VV_F)
	N = length(K)

	z = [zeros(2)]
	v = Array(Vector{Float64}, 0)
	for n = 1:N
		push!(v, -K[n]*z[n] - C[n])
		push!(z, A*z[n] + B*v[n])
	end
	return z, v
end
