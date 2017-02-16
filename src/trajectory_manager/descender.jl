######################################################################
# descender.jl
# provides descent engine
######################################################################
export Descender

export InverseStep
export InverseRootStep
export ArmijoLineSearch

type InverseStep <: Descender
	alpha::Float64
end

get_step_size(is::InverseStep, i::Int) = ir.alpha / i


type InverseRootStep <: Descender
	alpha::Float64
end

get_step_size(ir::InverseRootStep, i::Int) = ir.alpha / sqrt(i)


include("armijo_line_search.jl")
