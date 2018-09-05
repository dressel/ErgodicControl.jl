
export rrt_trajectory

struct Node
    depth::Int
    score::Float64
    parent::Int     # index of parent
    pose::Vector{Float64}
end

struct Graph
    open_nodes::Vector{Node}
    end_nodes::Vector{Node}
end
function Graph(x0::Vector{Float64})
    n_init = Node(1, 0.0, 0, x0)
    return Graph([n_init], Node[])
end


# K is the number of vertices
function rrt_trajectory(em::ErgodicManager, tm::TrajectoryManager; verbose::Bool=true, logging::Bool=false, K::Int=100, max_depth::Int=100, step_length::Float64=.1)
    # first create the graph
    G = Graph(tm.x0)

    best_depth = 1  # log the deepest we get (mostly for debugging)

    for k = 1:K
        # sample random point
        rx = rand() * em.domain.lengths[1] + em.domain.mins[1]
        ry = rand() * em.domain.lengths[2] + em.domain.mins[2]

        # find nearest node
        nearest_node = G.open_nodes[1]
        nearest_idx = 1
        idx = 1
        d2min = Inf
        for node in G.open_nodes
            dx = node.pose[1] - rx
            dy = node.pose[2] - ry
            d2 = dx*dx + dy*dy
            if d2 < d2min
                d2min = d2
                nearest_node = node
                nearest_idx = idx
            end
            idx += 1
        end

        # step from nearest node to new point
        dx = (rx - nearest_node.pose[1]) * step_length / sqrt(d2min)
        dy = (ry - nearest_node.pose[2]) * step_length / sqrt(d2min)
        new_pose = [nearest_node.pose[1] + dx, nearest_node.pose[2] + dy]

        # compute ergodic score
        traj = [new_pose, nearest_node.pose]
        current_node = nearest_node
        while current_node.parent > 0
            current_node = G.open_nodes[current_node.parent]
            push!(traj, current_node.pose)
        end
        es = ergodic_score(em, traj)

        # create new node and push into graph
        depth = nearest_node.depth + 1
        best_depth = max(depth, best_depth)
        new_node = Node(depth, es, nearest_idx, new_pose)
        if depth < max_depth
            push!(G.open_nodes, new_node)
        else
            push!(G.end_nodes, new_node)
        end

        if verbose
            println("num nodes = ", length(G.open_nodes))
            println("num closed = ", length(G.end_nodes))
        end
    end

    println("best_depth = ", best_depth)

    # now loop through the end nodes, finding the best one
    best_score = Inf
    best_node = G.open_nodes[1]     # just any node
    for node in G.end_nodes
        if node.score < best_score
            best_score = node.score
            best_node = node
        end
    end

    # make the full trajectory
    traj = [best_node.pose]
    current_node = best_node
    while current_node.parent > 0
        current_node = G.open_nodes[current_node.parent]
        push!(traj, current_node.pose)
    end

    # actually, reverse this
    return reverse(traj)
end
