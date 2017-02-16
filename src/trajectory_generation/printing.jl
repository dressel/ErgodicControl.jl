######################################################################
# printing.jl
######################################################################

function print_header()
	print_dashes()
	println(" iter  |ergodic score |control score |total score |direc deriv |step size")
	print_dashes()
end
function print_dashes()
	println("--------------------------------------------------------------------------")
end

function step_report(i::Int, es::Float64, cs::Float64, ts::Float64, dd::Float64, step_size::Float64)
	@printf " %-7i" i
	@printf " %-14.7f" es
	@printf " %-14.7f" cs
	@printf " %-12.7f" ts
	@printf " %-12.7f" dd
	@printf " %-7.5f" step_size
	println()
end
