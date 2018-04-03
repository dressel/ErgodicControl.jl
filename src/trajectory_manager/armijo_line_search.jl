######################################################################
# armijo_descender.jl
######################################################################
"""
type `ArmijoLineSearch <: Descender`

`ArmijoLineSearch(initial_step, c)`

Defaults are:
* `initial_step` = 10
* `c` = 0.5
* `max_iters` = 50
"""
mutable struct ArmijoLineSearch <: Descender
    initial_step::Float64
    c::Float64		# just a constant between 0 and 1
    max_iters::Float64

    function ArmijoLineSearch(initial_step::Real, c::Real, mi::Real)
        return new(float(initial_step), float(c), float(mi))
    end
    function ArmijoLineSearch(initial_step::Real, c::Real)
        return new(float(initial_step), float(c), 50.)
    end

    ArmijoLineSearch() = ArmijoLineSearch(10, 0.5, 50.)
end

function get_step_size(als::ArmijoLineSearch, em::ErgodicManager, tm::TrajectoryManager, xd::VVF, ud::VVF, zd::VVF, vd::VVF, ad::MF, bd::MF, K::Vector{MF}, i::Int)
    tau = 0.5
    step_size = als.initial_step

    # compute m = p' * grad f(x)
    m = directional_derivative(ad, bd, zd, vd)

    f_x = total_score(em, tm, xd, ud)

    xdn, udn = project(em, tm, K, xd, ud, zd, vd, step_size)
    ts = total_score(em, tm, xdn, udn)
    armijo_index = 0
    while (total_score(em, tm, xdn, udn) > f_x + step_size*als.c*m) && (armijo_index < als.max_iters)
        ts = total_score(em, tm, xdn, udn)
        step_size *= tau
        xdn, udn = project(em, tm, K, xd, ud, zd, vd, step_size)
        armijo_index += 1
    end
    return step_size
end


# kind of hacky, just for linear dynamics (no projection)
function get_step_size2(als::ArmijoLineSearch, em::ErgodicManager, tm::TrajectoryManager, xd::VVF, ud::VVF, zd::VVF, vd::VVF, ad::MF, bd::MF, K::Vector{MF}, i::Int)
    tau = 0.5
    step_size = als.initial_step

    # compute m = p' * grad f(x)
    m = directional_derivative(ad, bd, zd, vd)

    f_x = total_score(em, tm, xd, ud)

    xdn, udn = project2(em, tm, K, xd, ud, zd, vd, step_size)
    ts = total_score(em, tm, xdn, udn)
    armijo_index = 0
    while (total_score(em, tm, xdn, udn) > f_x + step_size*als.c*m) && (armijo_index < als.max_iters)
        ts = total_score(em, tm, xdn, udn)
        step_size *= tau
        xdn, udn = project2(em, tm, K, xd, ud, zd, vd, step_size)
        armijo_index += 1
    end
    return step_size
end
