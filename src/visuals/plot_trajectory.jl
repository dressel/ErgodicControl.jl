######################################################################
# plot_trajectory.jl
#
# Contains files for plotting the trajectory
######################################################################

# what other stuff do we need here?
# only marks, colors, etc
function plot_trajectory(xd::VVF;
                         ls::String="-",
                         lw::Real=1.0,
                         mfc::String="w",
                         ms::Real=9,
                         dims::Int=2
                        )

    N = length(xd)
    xvals = zeros(N)
    yvals = zeros(N)
    zvals = zeros(N)
    for i = 1:N
        xvals[i] = xd[i][1]
        yvals[i] = xd[i][2]
    end

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

    m = "."
    plot(xvals, yvals, "b", linestyle=ls, marker=m, lw=lw, ms=ms, mfc=mfc)

end

# plotting for multiple agents
# vtm is a vector of trajectory managers
function plot_trajectory(xd::VVF, vtm::VTM;
                         lw::Real=1,
                         ms::Real=9,
                         dims::Int=2
                        )

    num_agents = length(vtm)
    xds = vvf2vvvf(xd, vtm)
    for j = 1:num_agents
        plot_trajectory(xds[j], lw=lw, ms=ms,onlyMarks=onlyMarks,dims=dims)
    end
end
