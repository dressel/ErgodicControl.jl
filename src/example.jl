######################################################################
# example.jl
#
# includes some example problems
######################################################################

"""
`example(example_name::ASCIIString; K::Int=5, bins::Int=100)`

Possible strings are:
* single gaussian
* double gaussian
"""
function example(example_name::ASCIIString; K::Int=5, bins::Int=100)
	L = 1.0
	em = ErgodicManager(L, K, bins)

	if example_name == "single gaussian"
		mu = [L/2.0, L/2.0]
		Sigma = 0.03 * eye(2)
		phik!(em, mu, Sigma)
	elseif example_name == "double gaussian"
		# Create Gaussian distribution and its coefficients
		mu1 = [0.3, 0.7]
		Sigma1 = 0.025* eye(2)
		mu2 = [0.7, 0.3]
		Sigma2 = 0.025* eye(2)
		phik!(em, mu1, Sigma1, mu2, Sigma2)
	else
		println("example name not recognized")
	end
	return em
end
