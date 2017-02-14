######################################################################
# lqr.jl
#
# Solves LQR problems
######################################################################

export LQR
export apply_LQR_gains

function LQR(A::MF, B::MF, Q::MF, R::MF, N::Int)
	P = Array(Matrix{Float64}, N+1)
	G = Array(Matrix{Float64}, N)
	K = Array(Matrix{Float64}, N)
	P[N+1] = Q

	for n = (N-1):-1:0
		G[n+1] = R + (B' * P[n+1+1] * B)
		K[n+1] = inv(G[n+1]) * B' * P[n+1+1] * A
		P[n+1] = Q + (A' * P[n+1+1] * A) - (K[n+1]' * G[n+1] * K[n+1])
	end

	return K
end


function apply_LQR_gains(A::MF, B::MF, K::Vector{MF}, x0::Vector{Float64})
	x = [x0]
	u = Array(Vector{Float64}, 0)

	N = length(K)

	for n = 1:N
		push!(u, -K[n]*x[n])
		push!(x, A*x[n] + B*u[n])
	end
	return x, u
end
