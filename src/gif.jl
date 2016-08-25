######################################################################
# gif.jl
#
# makes gifs of ergodic trajectory generation
######################################################################

using Reel
export gif

"""
`gif(em::ErgodicManager, trajectories::Matrix{Float64}, num_trajectories::Int; show_score=true, fps::Int=20, right::Bool=true)`

Below is an example of how you might use this:

`xd, ud = clerc_trajectory(em, tm, max_iters=100, logging=true)`

`trajectories = readcsv("temp.csv")`

`gif(em, trajectories, 101, fps=17)`

The gif will be saved at `temp.gif`.
"""
function gif(em::ErgodicManager, trajectories::Matrix{Float64}, num_trajectories::Int; show_score=true, fps::Int=20, right::Bool=true)
	frames = Frames(MIME("image/png"), fps=fps)

	# ok, loop through all xd
	num_rows, n = size(trajectories)
	N = round(Int, num_rows / num_trajectories)

	# start with the first trajectory...
	xd = mat2traj(trajectories[1:N, :])
	plot(em, xd, show_score=show_score, right=right)
	push!(frames, gcf())
	close()

	# then do the rest
	for traj_idx = 2:num_trajectories
		start_idx = (traj_idx-1)*N + 1
		end_idx = start_idx + N - 1
		xd = mat2traj(trajectories[start_idx:end_idx,:])
		plot(em, xd, show_score=true)
		push!(frames, gcf())
		close()
	end

	write("temp.gif", frames)
end
