#= 
@authors: davidsantiagoquevedo, 17thSaint
Adapted from: https://github.com/17thSaint/finance-thesis/blob/master/Codes/fracworking-lang.jl

This script solves using Monte Carlo method the Fractional Langevin Equation:

d²x/dt² + dx/dt + D^alpha_t x = xi_white + xi_coloured 

* First guess:
Analytical solution of: d²x/dt² + D^alpha_t x = xi,
where, xi = xi_white + xi_coloured

x = x0 + v0*t + convolution(xi*t)

* Left-hand side for the residuals:
g(t) = d²x/dt² + dx/dt - xi = 0

* Right-hand side for the residuals:
D^alpha_t x

* Residual function:
d²x/dt² + dx/dt - xi + D^alpha_t x = 0
=#

DATA_PATH = "data/two_baths/"

include("../src/utils.jl")
include("../src/fBm_integration.jl")
include("../src/fractional_derivative.jl")
using .utils, .fbm_integration, .fractional_derivative, Plots, PyCall

h = parse(Float64,ARGS[1])
N = parse(Int64,ARGS[2])
T = parse(Int64,ARGS[3])
noise_id = parse(Int64,ARGS[4])

"""
    get_first_guess(times, noise, x0, v0)

First guess of the Monte Carlo Algorithm. 
Using the analytical solution of the fractional SDE: d²x/dt²  + D^alpha_t x = noise 
====> x = x0 + v0*t + convolution(noise*t)

ARGS
    times (Vector{Float64}): length of the configuration - number of points in the trajectory
	noise (Vector{Float64}): time delta between points
	x0 (float): initial position of the system
	v0 (float): initial velocity of the system
"""
function get_first_guess(times, noise, order, v0)
	py"""
	import numpy as np
	import sys
	sys.path.append("src/")
	import integration as itg
	import mittag_leffler as ml
	def solution_fle(t, noise, order, v0):
		t__ = np.array(t)
		noise__ = np.array(noise)
		def nonlinear_term(t):
			z = -t**(2-order)
			return t * ml.mittag_leffler_vector(z, 2-order, 2)
		conv = itg.convolution(nonlinear_term, noise__, t__)
		nonlinear = v0 * nonlinear_term(t__)
		return nonlinear + conv
	"""
    sol = py"solution_fle"(times, noise, order, v0)
    return sol
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
	lambda = 1
	beta = 1
	a = 1
	sliced_noise = noise[1:N-2]
    velocity = [(config[i+1] - config[i])/delta_t for i in 1:N-1]
	accel = [(velocity[i+1] - velocity[i])/delta_t for i in 1:N-2]
	
    g_of_t = lambda*accel + beta*lambda*velocity[1:N-2] - sliced_noise

	return g_of_t, velocity, accel
end

"""
    residuals(config, delta_t, noise, t_final, fd_order)
Get the residuals between the left-hand side and the right-hand side of the SDE.
ARGS
	config (array): current state/configuration of the system. For the langevin equation is the position
	delta_t (float): delta time between tow points on the configuration (N/t_fin)
	noise (array): noise term of the equation
    t_final (float): final time where the configuration is evaluated
    fd_order (float): order of the fractional derivative
"""
function residuals(config, delta_t, noise, t_final, fd_order)
	g_stuff = get_goft(config, delta_t, noise)
    N = length(config)
    frac_dev = grunwald_letnikov(fd_order, config, delta_t)[1:N-2]

	return abs.(g_stuff[1] + frac_dev)
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
    acc_rej_move(config, delta_t, noise, t_final, fd_order, chosen, step_size, metro_tol)
Accept or reject new configuration based on tolerance value (metro_tol)
ARGS
    config (array): current state/configuration of the system. For the langevin equation is the position
	delta_t (float): delta time between tow points on the configuration (N/t_fin)
	noise (array): noise term of the equation
    t_final (float): final time where the configuration is evaluated
    fd_order (float): order of the fractional derivative
    chosen (int): point in the configuration to shift
	step_size (array): maximum value to shift
    metro_tol (float): maximum tolerance value to accept or reject a shift
"""
function acc_rej_move(config, delta_t, noise, t_final, fd_order, chosen, step_size, metro_tol)
    N = length(config)
	start_resids = residuals(config, delta_t, noise, t_final, fd_order)
    shift_matrix = move_position(N, chosen, step_size)
	new_resids = residuals(config + shift_matrix, delta_t, noise, t_final, fd_order)
    exp_diff = exp.(new_resids - start_resids)
	
    checking = exp_diff[chosen-2] <= metro_tol # The acceptance rule must be evaluated in the position i-2 
												# because we are adjusting the second order derivative
    if checking
		return config + shift_matrix, 1
	else
		return config, 0
	end
end

function main_here(mc_steps, step_size, N, t_fin, noise, error_tol, metro_tol, first_guess, fd_order, sampling = false)
	delta_t = t_fin/N
	running_config = first_guess
	
	# Save configuration at every 10 attempted moves
	samp_freq = 10 
	samp_index = 1
	if sampling
		time_config = fill(0.0, (N, Int(mc_steps/samp_freq)))
		time_resids = fill(0.0, (N-2, Int(mc_steps/samp_freq)))
	else
		time_config = false
		time_resids = false
	end

	# Beginning of Monte Carlo loop
	num_wrong = [i for i in 3:N] # Assumed initial wrong points - all except for x0 and x1
	acc_rate = 0 # Acceptance rate of the Monte Carlo
	trial_index = 0

	for i in 1:mc_steps
		upper = 4
		if length(num_wrong) < 4
			upper = length(num_wrong)
		end
		
		for k in num_wrong[1:upper]
            #acc_rej_move(config, delta_t, noise, t_final, fd_order, chosen, step_size, metro_tol)
			movement = acc_rej_move(running_config, delta_t, noise, t_fin, fd_order, k, step_size, metro_tol)
			running_config = movement[1]
			acc_rate += movement[2]
			trial_index += 1
		end
			
		# Update wrong points
		num_wrong = []
		current_res = residuals(running_config, delta_t, noise, t_fin, fd_order)
		check_tol = [current_res[j] < error_tol for j in 1:N-2]
		for wr in 1:length(current_res)
			if ~ check_tol[wr]
				append!(num_wrong,[wr+2])
			end
		end

		# Save data every 10 steps
		if (i%samp_freq == 0) & (sampling)
			time_config[:, samp_index] = [running_config[x] for x in 1:N]
			time_resids[:, samp_index] = [current_res[y] for y in 1:N-2]
			samp_index += 1
		end

		# If every time point has residuals less than tolerance then solution is found
		if all(check_tol)
			println("Job $noise_id: Solution Found in $i Steps")
			return running_config, time_config, time_resids
		end
		
		# interface data
		if i%(mc_steps*0.01) == 0
			println("Job $noise_id:", " ", 100*i/mc_steps, "%, ", "Acceptance: ", acc_rate, "/", trial_index, ", Number Wrong: ", length(num_wrong))
			flush(stdout)
		end
	end
	
	println("Job $noise_id: No Solution")
	return running_config, time_config, time_resids
end

######################### ######################### #########################
#########################       RUN ALGORITHM       #########################
######################### ######################### #########################    

fd_order = 2 - 2*h
delta_t = T/N
noise_steps = 1
x0 = 0
v0 = 0

path_fbm = "$DATA_PATH"*"fBM-h-$h-$noise_id.hdf5"
path_obm = "$DATA_PATH"*"fBM-h-0.5-$noise_id.hdf5"
fbm = frac_brown_wiki2(h,N,T)
obm = frac_brown_wiki2(0.5,N,T)
write_data_hdf5(path_fbm, fbm)
write_data_hdf5(path_obm, obm)

w_noise = get_noise(0.5, noise_steps*N, T, noise_id, DATA_PATH)
c_noise = get_noise(h, noise_steps*N, T, noise_id, DATA_PATH)

noise = w_noise + c_noise
times = [i*T/N for i in 0:N-1]

fg = get_first_guess(times, noise, x0, v0)

error_tol = 0.0001
mc_steps = 1000000
metro_tol = 1.000001
step_size = 0.001#0.008*final_time/time_steps

println("Running job:$noise_id")
flush(stdout)
println("Warm-up run")
flush(stdout)
sol = main_here(10, step_size, N, T, noise, error_tol, metro_tol, fg, fd_order)
println("Sampling run")
flush(stdout)
sol = main_here(mc_steps, step_size, N, T, noise, error_tol, metro_tol, fg, fd_order)

write_data_hdf5(DATA_PATH * "two_bath-h-$h-noise$N.hdf5", (times, noise))
write_data_hdf5(DATA_PATH * "two_bath-h-$h-$N-v0$v0.hdf5", (times, sol[1]))