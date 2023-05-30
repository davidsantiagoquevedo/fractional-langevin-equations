module fractional_derivative
export grunwald_letnikov

using SpecialFunctions

function grunwald_letnikov(order, func, delta_t, factorial_top =  169)
	N = length(func)
	fract_deriv = [0.0 for _ in 0:N-1]
	fact_j = [exp(sum([log(i) for i in 1:j])) for j in 0:factorial_top]
	binom_part_j = gamma(order + 1)*[1/(fact_j[i]*gamma(order + 1 - (i-1))) for i in eachindex(fact_j)]
	
	for i in eachindex(fract_deriv) #ith point of the function
		sum_top = i
		if sum_top > factorial_top
			sum_top = factorial_top
		end
		for j in 0:sum_top-1 #discrete sum to calculate the integral
			func_here = func[i-j]
			binom_part = binom_part_j[j+1]
			fract_deriv[i] += ((-1)^(j))*binom_part*func_here/(delta_t^order)
		end
	end
	return fract_deriv
end
end