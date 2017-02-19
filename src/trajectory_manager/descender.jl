######################################################################
# descender.jl
# provides descent engine
######################################################################
export Descender

export InverseStep
export InverseRootStep
export ArmijoLineSearch
export ConstantStep

type InverseStep <: Descender
	alpha::Float64
end

function get_step_size(ir::InverseStep, em::ErgodicManager, tm::TrajectoryManager, xd::VV_F, ud::VV_F, zd::VV_F, vd::VV_F, ad::MF, bd::MF, K::Vector{MF}, i::Int)
	return ir.alpha / i
end


type InverseRootStep <: Descender
	alpha::Float64
end

#get_step_size(ir::InverseRootStep, i::Int) = ir.alpha / sqrt(i)
function get_step_size(ir::InverseRootStep, em::ErgodicManager, tm::TrajectoryManager, xd::VV_F, ud::VV_F, zd::VV_F, vd::VV_F, ad::MF, bd::MF, K::Vector{MF}, i::Int)
	return ir.alpha / sqrt(i)
end


type ConstantStep <: Descender
	alpha::Float64
end

function get_step_size(cs::ConstantStep, em::ErgodicManager, tm::TrajectoryManager, xd::VV_F, ud::VV_F, zd::VV_F, vd::VV_F, ad::MF, bd::MF, K::Vector{MF}, i::Int)
	return cs.alpha
end


include("armijo_line_search.jl")
