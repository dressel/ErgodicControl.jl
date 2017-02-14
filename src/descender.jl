######################################################################
# descender.jl
# provides descent engine
######################################################################
export Descender, InverseStep, InverseRootStep

abstract Descender

type InverseStep <: Descender
	alpha::Float64
end

get_step_size(is::InverseStep, i::Int) = ir.alpha / i


type InverseRootStep <: Descender
	alpha::Float64
end

get_step_size(ir::InverseRootStep, i::Int) = ir.alpha / sqrt(i)

