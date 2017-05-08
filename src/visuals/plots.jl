######################################################################
# plots.jl
######################################################################
using PyPlot
using PyCall
PyCall.PyDict(matplotlib["rcParams"])["font.family"]=["Times New Roman"]
PyCall.PyDict(matplotlib["rcParams"])["font.size"]=18

# must be called before plot
function axis_font(af::Real)
	PyCall.PyDict(matplotlib["rcParams"])["font.size"] = af
end

export plot, plot_trajectory, axis_font
"""
`plot(em::ErgodicManager, xd::VVF; alpha=1.0, cmap="Greys", show_score=true, right=true, lw=1.0, ms=6.0)`

`plot(em::ErgodicManager; alpha=1.0, cmap="Greys")`

The "Greys" cmap is dark where there is most density.
The "gray" cmap option is light where there is most density.

An `alpha` value closest to 1.0 is darker; less is more transparent.
"""
function plot(em::ErgodicManager, xd::VVF; alpha=1.0, cmap="Greys", show_score::Bool=true, right::Bool=true, lw::Float64=1.0, ms::Float64=6.0, onlyMarks::Bool=false, no_domain::Bool=false)
	plot_trajectory(xd, lw=lw, ms=ms, onlyMarks=onlyMarks)
	hold(true)
	plot(em, alpha=alpha, cmap=cmap, no_domain=no_domain)
	if show_score
		start_idx = right ? 1 : 0
		es = ergodic_score(em, xd, start_idx)
		title_string = "es = $(round(es,5))"
		title(title_string)
	end
end

function plot(em::ErgodicManagerR2; alpha=1.0, cmap="Greys",no_domain=false)
	a = [x_min(em), x_max(em), y_min(em), y_max(em)]
	if no_domain
		rar = zeros(size(em.phi'))
		imshow(rar, interpolation="none",cmap=cmap,origin="lower",extent=a,vmin=0, alpha=alpha)
	else
		imshow(em.phi', interpolation="none",cmap=cmap,origin="lower",extent=a,vmin=0, alpha=alpha)
	end
	axis(a)
end

function plot(em::ErgodicManagerSE2; alpha=1.0, cmap="Greys", no_domain=false)
	a = [x_min(em), x_max(em), y_min(em), y_max(em)]
	rar = size(em.phi)
	temp_phi = reshape(sum(em.phi,3), x_cells(em), y_cells(em))
	imshow(temp_phi', interpolation="none",cmap=cmap,origin="lower",extent=a,vmin=0, alpha=alpha)
	axis(a)
end

# assumes L = 1.0
# TODO: do I even use this anymore?
#function plot(mat::Matrix{Float64}; alpha=1.0, cmap="Greys")
#	L = 1.0
#	a = [0,L,0,L]
#	imshow(mat', interpolation="none",cmap=cmap,origin="lower",extent=a,vmin=0, alpha=alpha)
#	#labels()	# from an old package
#	axis(a)
#end

function plot(d::Domain, phi::Matrix{Float64}; alpha=1.0, cmap="Greys")
	a = [x_min(d), x_max(d), y_min(d), y_max(d)]
	imshow(phi', interpolation="none",cmap=cmap,origin="lower",extent=a,vmin=0, alpha=alpha)
	axis(a)
end

# what other stuff do we need here?
# only marks, colors, etc
function plot_trajectory(xd::VVF; lw::Real=1.0, ms::Real=6, onlyMarks=false)
	N = length(xd)
	xvals = zeros(N)
	yvals = zeros(N)
	for i = 1:N
		xvals[i] = xd[i][1]
		yvals[i] = xd[i][2]
	end
	ls = onlyMarks? "None" : "-"
	m = onlyMarks ? "o" : "."
	#PyPlot.plot(xvals, yvals, ".-", lw=lw, ms=ms)
	if onlyMarks
		#PyPlot.plot(xvals, yvals, linestyle=ls, marker=m, lw=lw, ms=ms, mfc="none")
		PyPlot.plot(xvals, yvals, linestyle=ls, marker=m, lw=lw, ms=ms, alpha=.1,mfc="black")
	else
		PyPlot.plot(xvals, yvals, linestyle=ls, marker=m, lw=lw, ms=ms)
	end
	#PyPlot.plot(xvals, yvals, line_style, lw=lw, ms=ms)
end

