######################################################################
# armijo_descender.jl
######################################################################
type ArmijoLineSearch <: Descender
	initial_step::Float64
	c::Float64		# just a constant between 0 and 1
	max_iters::Float64

	function ArmijoLineSearch(initial_step::Real, c::Real, mi::Real)
		return new(float(initial_step), float(c), float(mi))
	end
	function ArmijoLineSearch(initial_step::Real, c::Real)
		return new(float(initial_step), float(c), 50.)
	end

	ArmijoLineSearch() = ArmijoLineSearch(10, 0.1, 50.)

end

function get_step_size(als::ArmijoLineSearch, em::ErgodicManager, tm::TrajectoryManager, xd::VVF, ud::VVF, zd::VVF, vd::VVF, ad::MF, bd::MF, K::Vector{MF}, i::Int)
	tau = 0.5
	step_size = als.initial_step

	# compute m = p' * grad f(x)
	m = directional_derivative(ad, bd, zd, vd)

	f_x = total_score(em, tm, xd, ud)

	xdn, udn = project(em, tm, K, xd, ud, zd, vd, step_size)
	ts = total_score(em, tm, xdn, udn)
	#ts_arr = Float64[]
	#step_arr = Float64[]
	armijo_index = 0.0
	while (total_score(em, tm, xdn, udn) > f_x + step_size*als.c*m) && (armijo_index < als.max_iters)
		ts = total_score(em, tm, xdn, udn)
		#push!(ts_arr, ts)
		#push!(step_arr, step_size)
		step_size *= tau
		xdn, udn = project(em, tm, K, xd, ud, zd, vd, step_size)
		armijo_index += 1
	end
	#writecsv("ts.csv", ts_arr)
	#writecsv("step.csv", step_arr)
	return step_size
end
