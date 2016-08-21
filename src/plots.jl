######################################################################
# plots.jl
######################################################################
using PyPlot

export plot
function plot(em::ErgodicManager, xd::VV_F)
	N = length(xd)
	xvals = zeros(N)
	yvals = zeros(N)
	for i = 1:N
		xvals[i] = xd[i][1]
		yvals[i] = xd[i][2]
	end
	PyPlot.plot(xvals, yvals)
end
