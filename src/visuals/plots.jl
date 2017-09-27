######################################################################
# plots.jl
######################################################################
using PyPlot
import PyPlot.plot

# ensures the fonts will be Type 1 if we save as pdf
# if not, fonts might be Type 3, which most conferences don't allow
rc("text", usetex=true)

# specifying font family and size desires
rc("font", family="Times New Roman", size=18)

# must be called before plot
function axis_font(af::Real)
    rc("font", size=af)
end

# include helper files
include("plot_trajectory.jl")

# export my functions
export plot, plot_trajectory, axis_font

# export PyPlot functions user might want
export figure, savefig, xlabel, ylabel, zlabel, title, hold
export xlim, ylim, zlim


"""
`plot(em::ErgodicManager, xd::VVF; alpha=1.0, cmap="Greys", show_score=true, lw=1.0, ms=9)`

`plot(em::ErgodicManager; alpha=1.0, cmap="Greys")`

The "Greys" cmap is dark where there is most density.
The "gray" cmap option is light where there is most density.

An `alpha` value closest to 1.0 is darker; less is more transparent.
"""
function plot(em::ErgodicManager, xd::VVF;
              alpha=1.0,
              cmap="Greys",
              show_score::Bool=true,
              ls::String="-",
              lw::Real=1.5,
              mew::Real=1,
              mfc::String="w",
              ms::Real=10,
              no_domain::Bool=false
             )

    # If it is in R3, let the trajectory know
    dims = 2
    if typeof(em) == ErgodicManagerR3
        dims = 3
    end

    # plot the trajectory
    plot_trajectory(xd, ls="-", lw=lw, mew=mew, mfc=mfc, ms=ms, dims=dims)

    # plot domain
    plot(em, alpha=alpha, cmap=cmap, no_domain=no_domain)

    # determines if ergodic score should be shown
    if show_score
        es = ergodic_score(em, xd)
        title_string = "es = $(round(es,5))"
        title(title_string)
    end
end

# multi-agent version
# vtm is a vector of trajectory managers
function plot(em::ErgodicManager, xd::VVF, vtm::VTM;
              alpha=1.0,
              cmap="Greys",
              show_score::Bool=true,
              lw::Float64=1.0,
              ms::Float64=6.0,
              no_domain::Bool=false
             )

    xds = vvf2vvvf(xd, vtm)
    num_agents = length(vtm)

    # plot the trajectory
    for j = 1:num_agents
        xd = xds[j]
        plot_trajectory(xd, lw=lw, ms=ms)
    end

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
function plot(em::ErgodicManager, xd::VVF, n::Int;
              alpha=1.0,
              cmap="Greys",
              show_score::Bool=true,
              lw::Float64=1.0,
              ms::Float64=6.0,
              no_domain::Bool=false
             )

    # plotting the trajectory
    dims = 2
    plot_trajectory(xd[1:n], lw=lw, ms=ms, dims=dims)

    # plotting the distribution
    a = [x_min(em), x_max(em), y_min(em), y_max(em)]
    imshow(em.phi', interpolation="none",cmap=cmap,origin="lower",extent=a,vmin=0, alpha=alpha)
end

function plot(em::ErgodicManagerR2T, xd::VVF, n::Int;
              alpha=1.0,
              cmap="Greys",
              show_score::Bool=true,
              lw::Float64=1.0,
              ms::Float64=6.0,
              no_domain::Bool=false
             )

    # plotting the trajectory
    dims = 2
    plot_trajectory(xd[1:n], lw=lw, ms=ms, dims=dims)

    # plotting the distribution
    a = [x_min(em), x_max(em), y_min(em), y_max(em)]
    imshow(em.phi[:,:,n]', interpolation="none",cmap=cmap,origin="lower",extent=a,vmin=0, alpha=alpha)
end
function plot(em::ErgodicManagerR2T, xd::VVF, n::Int, vtm::VTM;
              alpha=1.0,
              cmap="Greys",
              show_score::Bool=true,
              lw::Float64=1.0,
              ms::Float64=6.0,
              no_domain::Bool=false
             )

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
    tick_params(direction="in")
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
    w = Weights(vec(em.phi))
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

function plot(em::ErgodicManagerSE2;
              alpha=1.0,
              cmap="Greys",
              no_domain=false
             )

    a = [x_min(em), x_max(em), y_min(em), y_max(em)]
    rar = size(em.phi)
    temp_phi = reshape(sum(em.phi,3), x_cells(em), y_cells(em))
    imshow(temp_phi', interpolation="none",cmap=cmap,origin="lower",extent=a,vmin=0, alpha=alpha)
    axis(a)
end


function plot(d::Domain, phi::Matrix{Float64}; alpha=1.0, cmap="Greys")
    a = [x_min(d), x_max(d), y_min(d), y_max(d)]
    imshow(phi', interpolation="none",cmap=cmap,origin="lower",extent=a,vmin=0, alpha=alpha)
    axis(a)
    tick_params(direction="in")
end
