#= 
@authors: davidsantiagoquevedo, 17thSaint
Adapted from: https://github.com/17thSaint/finance-thesis/blob/master/Codes/fracworking-lang.jl
=#

DATA_PATH = "data/"

include("../src/utils.jl")
include("../src/fBm_integration.jl")
include("../src/fractional_derivative.jl")
using .utils, .fbm_integration, .fractional_derivative, Plots, PyCall

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

"""
    get_goft(h, config, delta_t, noise, noise_steps, lambda, gam)
Integer order part of the differential equation

ARGS
	config (array): current state/configuration of the system. For the langevin equation is the position
	delta_t (float): delta time between tow points on the configuration (N/t_fin)
	noise (array): noise term of the equation
"""
function get_goft(config, delta_t, noise)
    N = length(config)

	sliced_noise = noise[1:N-2]
    velocity = [(config[i+1] - config[i])/delta_t for i in 1:N-1]
	accel = [(velocity[i+1] - velocity[i])/delta_t for i in 1:N-1]
	
    g_of_t = accel - sliced_noise

	return g_of_t, velocity, accel
end

"""
	get_resids(config, delta_t, noise)
Get the residuals between the left-hand side and the right-hand side of the SDE.
ARGS
	config (array): current state/configuration of the system. For the langevin equation is the position
	delta_t (float): delta time between tow points on the configuration (N/t_fin)
	noise (array): noise term of the equation
    fd_order (float): order of the fractional derivative
"""
function residuals(config, delta_t, noise, fd_order)
	g_stuff = get_goft(config, delta_t, noise)
    N = length(config)
    frac_dev = grunwald_letnikov(fd_order, N, delta_t, config)[1:N-2]

	return abs.(g_stuff[1] - frac_dev)
end

"""
	move_position(N, chosen, step_size)
Suggest a shift matrix to change the configuration at a certain point in time using a random number with a maximum step size
ARGS
	N (int): length of the configuration - number of points in the trajectory
	chosen (int): point in the configuration to shift
	step_size (array): maximum value to shift
"""
function move_position(N, chosen, step_size)
	shift_matrix = [0.0 for _ in 1:N]
	shift_matrix[chosen] += rand(-1:2:1)*rand(Float64)*step_size
	return shift_matrix
end

"""
	acc_rej_move(config, chosen, step_size, delta_t, noise, fd_order, metro_tol)
Accept or reject new configuration based on tolerance value (metro_tol)
ARGS
    config (array): current state/configuration of the system. For the langevin equation is the position
	chosen (int): point in the configuration to shift
	step_size (array): maximum value to shift
	delta_t (float): delta time between tow points on the configuration (N/t_fin)
	noise (array): noise term of the equation
	metro_tol (float): maximum tolerance value to accept or reject a shift
    fd_order (float): order of the fractional derivative
"""
function acc_rej_move(config, chosen, step_size, delta_t, noise, fd_order, metro_tol)
    N = lenght(config)
	start_resids = residuals(config, delta_t, noise, fd_order)
    shift_matrix = move_position(N, chosen, step_size)
	new_resids = residuals(config + shift_matrix, delta_t, noise, fd_order)
    exp_diff = exp.(new_resids - start_resids)
	
    checking = exp_diff[chosen-2] <= metro_tol # The acceptance rule must be evaluated in the position i-2 
												# because we are adjusting the second order derivative
    if checking
		return config + shift_matrix, 1
	else
		return config, 0
	end
end

######################### ######################### #########################
#########################       RUN ALGORITHM       #########################
######################### ######################### #########################    

h = 0.5
final_time = 30
N = 500
noise_steps = 1
x0 = 0
v0 = 0
noise = get_noise(h, noise_steps*N, final_time, 1, DATA_PATH)
times = [i*final_time/N for i in 0:N-1]

fg = get_first_guess(times, noise, x0, v0)

