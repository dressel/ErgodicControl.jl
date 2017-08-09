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

# Export my functions and export PyPlot functions
export plot, plot_trajectory, axis_font
export figure, savefig, xlabel, ylabel, zlabel, title, hold
export xlim, ylim, zlim
"""
`plot(em::ErgodicManager, xd::VVF; alpha=1.0, cmap="Greys", show_score=true, lw=1.0, ms=6.0)`

`plot(em::ErgodicManager; alpha=1.0, cmap="Greys")`

The "Greys" cmap is dark where there is most density.
The "gray" cmap option is light where there is most density.

An `alpha` value closest to 1.0 is darker; less is more transparent.
"""
function plot(em::ErgodicManager, xd::VVF; alpha=1.0, cmap="Greys", show_score::Bool=true, lw::Float64=1.0, ms::Float64=6.0, onlyMarks::Bool=false, no_domain::Bool=false)

	# If it is in R3, let the trajectory know
	dims = 2
	if typeof(em) == ErgodicManagerR3
		dims = 3
	end

	# plot the trajectory
	plot_trajectory(xd, lw=lw, ms=ms, onlyMarks=onlyMarks, dims=dims)

	# hold and plot domain
	hold(true)
	plot(em, alpha=alpha, cmap=cmap, no_domain=no_domain)

	# determines if ergodic score should be shown
	if show_score
		es = ergodic_score(em, xd)
		title_string = "es = $(round(es,5))"
		title(title_string)
	end
end

# multi-agent version
function plot(em::ErgodicManager, xd::VVF, vtm::Vector{TrajectoryManager}; alpha=1.0, cmap="Greys", show_score::Bool=true, lw::Float64=1.0, ms::Float64=6.0, onlyMarks::Bool=false, no_domain::Bool=false)

	xds = vvf2vvvf(xd, vtm)
	num_agents = length(vtm)

	# plot the trajectory
	for j = 1:num_agents
		xd = xds[j]
		plot_trajectory(xd, lw=lw, ms=ms, onlyMarks=onlyMarks)
	end

	# hold and plot domain
	hold(true)
	plot(em, alpha=alpha, cmap=cmap, no_domain=no_domain)

	# determines if ergodic score should be shown
	#if show_score
	#	es = ergodic_score(em, xd)
	#	title_string = "es = $(round(es,5))"
	#	title(title_string)
	#end
end

# special plotting function for R2T
# input n is a Julia index; it goes from 1 to N+1 (instead of 0 to N)
function plot(em::ErgodicManager, xd::VVF, n::Int; alpha=1.0, cmap="Greys", show_score::Bool=true, lw::Float64=1.0, ms::Float64=6.0, onlyMarks::Bool=false, no_domain::Bool=false)

	# plotting the trajectory
	dims = 2
	plot_trajectory(xd[1:n], lw=lw, ms=ms, onlyMarks=onlyMarks, dims=dims)

	# plotting the distribution
	a = [x_min(em), x_max(em), y_min(em), y_max(em)]
	imshow(em.phi', interpolation="none",cmap=cmap,origin="lower",extent=a,vmin=0, alpha=alpha)
end
function plot(em::ErgodicManagerR2T, xd::VVF, n::Int; alpha=1.0, cmap="Greys", show_score::Bool=true, lw::Float64=1.0, ms::Float64=6.0, onlyMarks::Bool=false, no_domain::Bool=false)

	# plotting the trajectory
	dims = 2
	plot_trajectory(xd[1:n], lw=lw, ms=ms, onlyMarks=onlyMarks, dims=dims)

	# plotting the distribution
	a = [x_min(em), x_max(em), y_min(em), y_max(em)]
	imshow(em.phi[:,:,n]', interpolation="none",cmap=cmap,origin="lower",extent=a,vmin=0, alpha=alpha)
end
function plot(em::ErgodicManagerR2T, xd::VVF, n::Int, vtm::VTM; alpha=1.0, cmap="Greys", show_score::Bool=true, lw::Float64=1.0, ms::Float64=6.0, onlyMarks::Bool=false, no_domain::Bool=false)

	# plotting the trajectory
	dims = 2
	plot_trajectory(xd[1:n], vtm, lw=lw, ms=ms, onlyMarks=onlyMarks, dims=dims)

	# plotting the distribution
	a = [x_min(em), x_max(em), y_min(em), y_max(em)]
	imshow(em.phi[:,:,n]', interpolation="none",cmap=cmap,origin="lower",extent=a,vmin=0, alpha=alpha)
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

# TODO: this is not ready to go yet
function plot(em::ErgodicManagerR3; alpha=1.0, cmap="Greys",no_domain=false)
	#if no_domain
	#	rar = zeros(size(em.phi'))
	#	imshow(rar, interpolation="none",cmap=cmap,origin="lower",extent=a,vmin=0, alpha=alpha)
	#else
	#	imshow(em.phi', interpolation="none",cmap=cmap,origin="lower",extent=a,vmin=0, alpha=alpha)
	#end

	# this is super sketch but let's try it
	x_temp, y_temp, z_temp = get_pts2(em)
	scatter3D(x_temp,y_temp,z_temp)
	xlim(x_min(em), x_max(em))
	ylim(y_min(em), y_max(em))
	zlim(z_min(em), z_max(em))
end

function get_pts2(em::ErgodicManagerR3)
	w = WeightVec(vec(em.phi))
	nx,ny,nz = x_cells(em), y_cells(em), z_cells(em)
	x_temp, y_temp, z_temp = Float64[], Float64[], Float64[]
	xs,ys,zs = x_size(em), y_size(em), z_size(em)

	for i = 1:1000
		xi, yi, zi = ind2sub(size(em.phi), sample(w))
		push!(x_temp, (xi-.5) * xs)
		push!(y_temp, (yi-.5) * ys)
		push!(z_temp, (zi-.5) * zs)
	end
	return x_temp, y_temp, z_temp
end

function get_pts(em::ErgodicManagerR3)
	x_temp = Float64[]
	y_temp = Float64[]
	z_temp = Float64[]
	for xi = 1:x_cells(em)
		for yi = 1:y_cells(em)
			for zi = 1:z_cells(em)
				if em.phi[xi,yi,zi] > 0.5
					push!(x_temp, (xi - .5)*x_size(em))
					push!(y_temp, (yi - .5)*y_size(em))
					push!(z_temp, (zi - .5)*z_size(em))
				end
			end
		end
	end
	return x_temp, y_temp, z_temp
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
function plot_trajectory(xd::VVF; lw::Real=1.0, ms::Real=6, onlyMarks=false, dims::Int=2)
	N = length(xd)
	xvals = zeros(N)
	yvals = zeros(N)
	zvals = zeros(N)
	for i = 1:N
		xvals[i] = xd[i][1]
		yvals[i] = xd[i][2]
	end
	ls = onlyMarks? "None" : "-"
	m = onlyMarks ? "o" : "."

	# if in 3 dimensions
	if dims == 3
		for i = 1:N
			zvals[i] = xd[i][3]
		end
		plot3D(xvals, yvals, zvals, ls=ls, marker=m, ms=ms, lw=lw)
		xlabel(L"x")
		ylabel(L"y")
		zlabel(L"z")
		return
	end

	#PyPlot.plot(xvals, yvals, ".-", lw=lw, ms=ms)
	if onlyMarks
		#PyPlot.plot(xvals, yvals, linestyle=ls, marker=m, lw=lw, ms=ms, mfc="none")
		PyPlot.plot(xvals, yvals, linestyle=ls, marker=m, lw=lw, ms=ms, alpha=.1,mfc="black")
	else
		PyPlot.plot(xvals, yvals, linestyle=ls, marker=m, lw=lw, ms=ms)
	end
	#PyPlot.plot(xvals, yvals, line_style, lw=lw, ms=ms)
end

function plot_trajectory(xd::VVF, vtm::VTM; lw::Real=1.0, ms::Real=6, onlyMarks=false, dims::Int=2)

	num_agents = length(vtm)
	xds = vvf2vvvf(xd, vtm)
	hold(true)
	for j = 1:num_agents
		plot_trajectory(xds[j], lw=lw, ms=ms,onlyMarks=onlyMarks,dims=dims)
	end
end
