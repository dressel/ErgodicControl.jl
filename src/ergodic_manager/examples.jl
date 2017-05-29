######################################################################
# examples.jl
# 
# This used to be in constructors of ergodic managers
# It added a lot of clutter and I didn't like it there
######################################################################
function ErgodicManagerR2(example_name::String; K::Int=5, bins::Int=100)
	L = 1.0
	d = Domain([0.,0.], [L, L], [bins,bins])

	if example_name == "single gaussian"
		mu = [L/2.0, L/2.0]
		Sigma = 0.03 * eye(2)
		phi = gaussian(d, mu, Sigma)
		return ErgodicManagerR2(d, phi, K)
	elseif example_name == "double gaussian"
		mu1 = [0.3, 0.7]
		Sigma1 = 0.025 * eye(2)
		mu2 = [0.7, 0.3]
		Sigma2 = 0.025 * eye(2)
		#phi = gaussian(d, [mu1,mu2], [Sigma1, Sigma2], [.5,.5])
		phi = gaussian(d, [mu1,mu2], [Sigma1, Sigma2])
		return ErgodicManagerR2(d, phi, K)
	else
		error("example name not recognized")
	end
end

function ErgodicManagerSE2(example_name::String; K::Int=5,bins::Int=50)
	L = 1.0
	d = Domain([0.,0.,-pi], [L,L,pi], [bins,bins,bins])
	println("h2")
	println("h3")

	if example_name == "single gaussian"
		mu = [L/2.0, L/2.0, 0]
		Sigma = 0.03 * eye(3)
		phi = gaussian(d, mu, Sigma)
		return ErgodicManagerSE2(d, phi, K)

	elseif example_name == "double gaussian"
		# How I should do it...
		mu1 = [0.3, 0.7]
		Sigma1 = 0.025* eye(2)
		mu2 = [0.7, 0.3]
		Sigma2 = 0.025* eye(2)
		phi = gaussian(d, mu1, Sigma1, mu2, Sigma2)
		return ErgodicManagerR2(d, phi, K)

		# Create Gaussian distribution and its coefficients
		#mu1 = [0.3, 0.7]
		#Sigma1 = 0.025* eye(2)
		#mu2 = [0.7, 0.3]
		#Sigma2 = 0.025* eye(2)
		#phi!(em, mu1, Sigma1, mu2, Sigma2)
		#decompose!(em)
	else
		error("example name not recognized")
	end
	return em
end
