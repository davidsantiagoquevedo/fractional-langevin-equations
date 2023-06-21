#= 
@authors: davidsantiagoquevedo, 17thSaint
Adapted from: https://github.com/17thSaint/finance-thesis/blob/master/Codes/fracworking-lang.jl
=#

DATA_PATH = "data/"
FIG_PATH = "analysis/inspect_benchmark2_fle_files/"

include("../src/utils.jl")
include("../src/fBm_integration.jl")
include("../src/fractional_derivative.jl")
using .utils, .fractional_derivative, Plots, PyCall

"""
    get_first_guess(times, noise, x0, v0)

First guess of the Monte Carlo Algorithm. 
Using the analytical solution of the integer SDE: d^2x/dt = noise 
====> x = x0 + v0*t + convolution(noise*t)

ARGS
    times (Vector{Float64}): length of the configuration - number of points in the trajectory
	noise (Vector{Float64}): time delta between points
	x0 (float): initial position of the system
	v0 (float): initial velocity of the system
"""
function get_first_guess(times, noise, x0, v0)
    py"""
    import sys
    sys.path.append("src/")
    import integration as itg
    import numpy as np

    def line(t):
        return t

    def solution_sde(t, noise, x0, v0):
        t_ = np.array(t)
        noise_ = np.array(noise)
        return x0 + v0*t_ + itg.convolution(line, noise_, t_)
    """
    a = py"solution_sde"(times, noise, x0, v0)
    return a
end

function analytical(times, noise, order, v0)
	py"""
	import numpy as np
	import sys
	sys.path.append("src/")
	import integration as itg
	import mittag_leffler as ml
	def solution_fle_white(t, noise, order, v0):
		t__ = np.array(t)
		noise__ = np.array(noise)
		def nonlinear_term(t):
			z = -t**(2-order)
			return t * ml.mittag_leffler_vector(z, 2-order, 2)
		conv = itg.convolution(nonlinear_term, noise__, t__)
		nonlinear = v0 * nonlinear_term(t__)
		return nonlinear + conv
	"""
    sol = py"solution_fle_white"(times, noise, order, v0)
    return sol
end

#######    

h = 0.4
fd_order = 2 - 2*h

final_time = 30
N = 350
delta_t = final_time/N
noise_steps = 1
x0 = 0
v0 = 0

percentage = 1

noise = read_hdf5(DATA_PATH * "fle-h-$h-noise$N.hdf5")[2]
times = read_hdf5(DATA_PATH * "fle-h-$h-$N-v0$v0.hdf5")[1]
sol = read_hdf5(DATA_PATH * "fle-h-$h-$N-v0$v0.hdf5")[2]

fg = get_first_guess(times, noise, x0, v0)
anl = analytical(times, noise, fd_order, fg[2]/delta_t)

plot(times[1:Int(N*percentage)], sol[1:Int(N*percentage)], label = "MC - Out", marker =:diamond)
plot!(times[1:Int(N*percentage)], anl[1:Int(N*percentage)], label = "Analytical")
#png("$FIG_PATH/fle_h=$h-N$N-anl") 