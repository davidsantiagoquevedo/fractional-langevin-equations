#=
@authors: 17thsaint, davidsantiagoquevedo
Adapted from: https://github.com/17thSaint/finance-thesis
=#

module fbm_integration
export frac_brown_wiki2, get_noise

include("utils.jl")
using .utils, Cubature, SpecialFunctions, HypergeometricFunctions

function frac_brown_wiki2(h, n, t_fin, print_satus = true)
    #=
    generator of fractional brownian trajectories based on Stochastic integral definition
    See Method 2: https://en.wikipedia.org/wiki/Fractional_Brownian_motion
	
    h: hurst number
    n: number of points of the path
    t_fin: final time (T) of the trajectory
    =#
	times = [i*t_fin/n for i in 0:n]
	dB = [rand(Float64)*rand(-1:2:1)*sqrt(t_fin/n) for i in 1:n]
	bh = fill(0.0,n)
	for j in 1:n
		if print_satus
			if j%(0.05*n) == 0
				println("H=",h,", ",100*j/n,"%",", fBM")
			end
		end 
		for i in 0:j-1
			function integrand(s)
				return (((times[j+1]-s)^(h-0.5))/gamma(h+0.5))*_₂F₁(h-0.5,0.5-h,h+0.5,1-times[j+1]/s)
			end
			part = hquadrature(integrand,times[i+1],times[i+2])[1]*dB[i+1]*n/t_fin
			bh[j] += part
		end
	end
	full_bh = append!([0.0],bh)
	return times, full_bh
end

function get_noise(h, N, t_fin, which = rand(1:10), dir ="", make_new = false)
	fBM = read_hdf5_data(h, which, dir, true, N+1)[2]
	if make_new
		fBM = frac_brown_wiki2(h, N, t_fin)[2]
	end

	return [t_fin*(fBM[i+1] - fBM[i])/N for i in 1:N]
end

end