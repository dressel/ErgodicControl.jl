######################################################################
# armijo_descender.jl
######################################################################
type ArmijoLineSearch <: Descender
end

function get_step_size(als::ArmijoLineSearch, em::ErgodicManager, tm::TrajectoryManager, xd::VV_F, ud::VV_F, zd::VV_F, vd::VV_F, ad::MF, bd::MF, K::Vector{MF}, i::Int)
	#print("armijo starting...")
	tau = 0.5
	c = 0.5
	#c = 0.999
	#c = 0.0001
	c = 0.5
	c = 1e-6
	step_size = 0.01

	# compute m = p' * grad f(x)
	#m = directional_derivative(ad, bd, zd, vd)
	#normalizer!(zd, vd)
	m = directional_derivative(ad, bd, zd, vd)
	#println("m = ", m)
	#m = scaled_dd(ad, bd, zd, vd)

	f_x = total_score(em, tm, xd, ud)
	#println("f_x = ", f_x)
	#println("m = ", m)

	xdn, udn = project(em, tm, K, xd, ud, zd, vd, step_size)
	ts = total_score(em, tm, xdn, udn)
	#println("ts = ", round(ts,3))
	ts_arr = Float64[]
	step_arr = Float64[]
	while total_score(em, tm, xdn, udn) > f_x + step_size*c*m
		#print("step_size = ", step_size)
		#print(", f_x + step_size*c*m = ", round(f_x + step_size*c*m,3))
		ts = total_score(em, tm, xdn, udn)
		push!(ts_arr, ts)
		push!(step_arr, step_size)
		#println(", ts = ", round(ts,3))
		step_size *= tau
		xdn, udn = project(em, tm, K, xd, ud, zd, vd, step_size)
	end
	writecsv("ts.csv", ts_arr)
	writecsv("step.csv", step_arr)
	#println("step_size = ", step_size)
	return step_size
end
