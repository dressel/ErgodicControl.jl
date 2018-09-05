######################################################################
# trajectory_generation.jl
######################################################################

# helper files
include("gradients.jl")
include("scoring.jl")
include("projection.jl")
include("printing.jl")

# actual methods
include("pto.jl")
include("rrt.jl")
include("pto_multi.jl")
include("pto_linear.jl")
include("smc.jl")
include("smc_multi.jl")
include("max.jl")
include("kmeans.jl")
#include("rig.jl")

# some additional helpers
include("info.jl")


function check_convergence(es::Float64, es_crit::Float64, i::Int, max_iters::Int, dd::Float64, dd_crit::Float64, verbose::Bool, es_count::Int)
    not_finished = true
    if es < es_crit
        not_finished = false
        if verbose
            println("reached ergodic criterion...")
        end
    end
    if i > max_iters
        not_finished = false
        if verbose
            println("max iterations reached...")
        end
    end
    if abs(dd) < dd_crit
        not_finished = false
        if verbose
            println("reached directional derivative criterion...")
        end
    end
    if es_count > 50
        not_finished = false
        if verbose
            println("We've been stuck for 50 iterations...")
        end
    end
    return not_finished
end

# called if logging, not meant for general use
function save(outfile::IOStream, xd::VVF)
    n = length(xd[1])
    for xi in xd
        for i = 1:(n-1)
            wi = xi[i]
            write(outfile,"$(xi[i]),")
        end
        write(outfile,"$(xi[n])\n")
    end
end
