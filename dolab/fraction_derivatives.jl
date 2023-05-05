function get_fractderiv(order, steps, delta_t, func)
    fract_deriv = [0.0 for _ in 1:steps-1]


	fact_j = exp(sum([log(i) for i in 1:j]))
	binom_part = gamma(order + 1)/(fact_j*gamma(order + 1 - j))

	for i in 1:steps-1
		top = i
		if i >= 169
			top = 169
		end
		for j in 0:top
			func_here = og_func[i-j+1]
			binom_part = binom_part1[j+1]
			fract_deriv[i] += ((-1)^(j))*binom_part*func_here/(delta_t/fracderiv_steps)^(2-2*h)
		end
	end
	return fract_deriv
end


time_told = 200
factorial_top = time_told
if time_told > 169
    factorial_top = 169
end
fact_j = [exp(sum([log(i) for i in 1:j])) for j in 0:factorial_top]
binom_part1 = gamma(3-2*h).*[1/(fact_j[i]*gamma(3-2*h-(i-1))) for i in 1:length(fact_j)]


#get_fractderiv(h,delta_t,steps,og_func,f0,noise_steps,binom_part1,cap=1)
#get_fractderiv(h,delta_t,steps,config,0.0,noise_steps,binom_part1)