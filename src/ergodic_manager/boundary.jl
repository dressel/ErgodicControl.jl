######################################################################
# boundary.jl
#######################################################################
export PointsBoundary
export in_boundary
export find_closest_boundary, interpolate_points


##########################################
# Points Boundary
##########################################

mutable struct PointsBoundary
	# user provided
	points::Vector{Tuple{Float64, Float64}} # [(x1,y1), (x2,y2)....]

	function PointsBoundary(points)
        b = new()
        b.points = points
		return b
	end
    
end

function is_point_inside_set(p, points::Vector{Tuple{Float64,Float64}}, tol::Float64=1e-6)
    # Count the number of intersections between the horizontal line passing through
    # point p and the polygon defined by the set of points.
    if points == Tuple{Float64, Float64}[]
        return true
    end

    num_intersections = 0
    for i in 1:length(points)
        j = i % length(points) + 1
        # Check if the line intersects the edge defined by points i and j.
        if (points[i][2] > p[2] + tol) != (points[j][2] > p[2] + tol)
            # Compute the x-coordinate of the intersection point using the
            # equation of the line passing through points i and j.
            t = (p[2] - points[i][2]) / (points[j][2] - points[i][2])
            x_int = t * (points[j][1] - points[i][1]) + points[i][1]
            # If the intersection point is to the right of the query point, count it as an intersection.
            if x_int > p[1] + tol
                num_intersections += 1
            end
        end
    end
    # Return true if the number of intersections is odd (i.e., the point is inside
    # the polygon) and false otherwise (i.e., the point is outside the polygon).
    return isodd(num_intersections)
end

function in_boundary(p, boundary::PointsBoundary)   
    return is_point_inside_set(p, boundary.points)
end

function find_closest_boundary(p, points::Vector{Tuple{Float64,Float64}}, tol::Float64=1e-6)
    # Initialize the minimum distance and closest point
    min_dist = Inf
    closest_point = nothing
    
    # Iterate over each point in the set
    for point in points
        # Compute the Euclidean distance between the query point and the current point
        dist = sqrt((p[1] - point[1])^2 + (p[2] - point[2])^2)
        
        # If the distance is smaller than the current minimum and greater than the tolerance, update the minimum distance and closest point
        if (dist < min_dist) && (dist > tol)
            min_dist = dist
            closest_point = point
        end
    end
    
    # Return the closest point, or nothing if no points were found within the tolerance
    return closest_point
end

function find_closest_boundary(p, boundary::PointsBoundary)
    return find_closest_boundary(p, boundary.points)
end

function sort_points(points::Vector{Tuple{Float64, Float64}})
    # Compute the centroid of the points
    cx, cy = sum([p[1] for p in points]) / length(points), sum([p[2] for p in points]) / length(points)
    # Define a function to compute the angle between two points and the x-axis
    angle(x::Tuple{Float64, Float64}) = atan(x[2] - cy, x[1] - cx)
    # Sort the points by angle with respect to the centroid
    sorted_points = sort(points, by=angle)

    return sorted_points
end

function remove_duplicates(points::Vector{Tuple{Float64,Float64}})
    seen = Vector{Tuple{Float64,Float64}}()
    unique_points = Vector{Tuple{Float64,Float64}}()
    for point in points
        if !(point in seen)
            push!(unique_points, point)
            push!(seen, point)
        else
            unique_points = filter(p -> p != point, unique_points)
        end
    end
    return unique_points
end

function interpolate_points(points::Vector{Tuple{Float64, Float64}}, N::Int)
    new_points = Vector{Tuple{Float64, Float64}}(undef, 0)

    for i in 1:length(points)-1
        start_point = points[i]
        end_point = points[i+1]
        step = (end_point .- start_point) ./ (N+1)

        for j in 0:N
            interpolated_point = start_point .+ j .* step
            push!(new_points, interpolated_point)
        end
    end

    # Include the last point from the original set
    push!(new_points, points[end])

    return new_points
end

# function convert_poly_boundary2_points_boundary(boundary::Boundary, poly_side)
#     points = Tuple{Float64,Float64}[]
#     for rect in boundary.rectangles
#         p1x = rect[1]
#         p1y = rect[2]
#         p2x = rect[3]
#         p2y = rect[4]

#         push!(points, (p1x, p1y))
#         push!(points, (p1x, p1y+poly_side))
#         push!(points, (p1x+poly_side, p1y+poly_side))
#         push!(points, (p1x+poly_side, p1y))
#     end

#     sorted_points = sort_points(remove_duplicates(points))
#     push!(sorted_points, sorted_points[1]) # connect the last point to the first point

#     interp_sorted_points = interpolate_points(sorted_points, N_boundary_interp)

#     return PointsBoundary(interp_sorted_points)
# end

