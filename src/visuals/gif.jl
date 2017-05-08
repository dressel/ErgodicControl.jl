######################################################################
# gif.jl
#
# makes gifs of ergodic trajectory generation
######################################################################

using Reel
export gif

"""
`gif(em::ErgodicManager, tm::TrajectoryManager, trajectory_file::String="temp.csv"; show_score=true, fps::Int=20, right::Bool=true)`

Below is an example of how you might use this. You `must` set `logging` to true during the trajectory generation.

`xd, ud = pto_trajectory(em, tm, logging=true)`

`gif(em, tm, fps=15)`

The gif will be saved at `temp.gif`.
"""
function gif(em::ErgodicManager, tm::TrajectoryManager, trajectory_file::String="temp.csv"; show_score=true, fps::Int=17, right::Bool=true)
	# First, let's read the trajectories in...
	trajectories = readcsv(trajectory_file)


	frames = Frames(MIME("image/png"), fps=fps)

	# ok, loop through all xd
	num_rows, n = size(trajectories)
	#N = round(Int, num_rows / num_trajectories)
	N = tm.N + 1
	num_trajectories = round(Int, num_rows / N)

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
