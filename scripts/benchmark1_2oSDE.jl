#= 
@authors: davidsantiagoquevedo, 17thSaint
Adapted from: https://github.com/17thSaint/finance-thesis/blob/master/Codes/fracworking-lang.jl
=#

DATA_PATH = "data/"

include("../src/utils.jl")
include("../src/fBm_integration.jl")
using .utils, .fbm_integration, PyPlot
pygui(true)

"""
    get_first_guess(N, delta_t, x0, v0)

First guess of the Monte Carlo Algorithm. 
For the Fractional Langevin equations, the first guess is the analytical solution of the integer order part of de differential equation.
However, for this benchmark we use x = [x0, x1, 0, 0, 0, ...]

ARGS
    N (int): length of the configuration - number of points in the trajectory
	delta_t (float): time delta between points
	x0 (float): initial position of the system
	v0 (float): initial velocity of the system
"""
function get_first_guess(N, delta_t, x0, v0)
	x1 = x0 + v0*delta_t
	x_guess = [0.0 for _ in 3:N]
	x_guess = append!([x0, x1], x_guess)

    return x_guess
end

"""
    get_goft(h, config, delta_t, noise, noise_steps, lambda, gam)
Integer order part of the differential equation
Removed get_scale_inv_vals: Notice that in the original implementation they were returning 1.

ARGS
	config (array): current state/configuration of the system. For the langevin equation is the position
	delta_t (float): delta time between tow points on the configuration (N/t_fin)
	noise (array): noise term of the equation
"""
function get_goft(config, delta_t, noise)
	sliced_noise = [noise[i] for i in 1:length(config)-2]
	velocity = [(config[i+1] - config[i])/delta_t for i in 1:length(config)-1]
	accel = [(velocity[i+1] - velocity[i])/delta_t for i in 1:length(velocity)-1]
	g_of_t = accel - sliced_noise

	return g_of_t, velocity, accel
end

"""
	get_resids(config, delta_t, noise)
Get the residuals between the left-hand side and the right-hand side of the SDE.
For the fractional Langevin equation it's the residual between the analytically tratable part and the fractional derivative.
ARGS
	config (array): current state/configuration of the system. For the langevin equation is the position
	delta_t (float): delta time between tow points on the configuration (N/t_fin)
	noise (array): noise term of the equation
"""
function get_resids(config, delta_t, noise)
	g_stuff = get_goft(config, delta_t, noise)
	right_hand_side = [0.0 for _ in 1:length(g_stuff[1])]

	return abs.(g_stuff[1] - right_hand_side)
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
	acc_rej_move(config, N, chosen, step_size, delta_t, noise, metro_tol)
Accept or reject new configuration based on tolerance value (metro_tol)
ARGS
	N (int): length of the configuration - number of points in the trajectory
	chosen (int): point in the configuration to shift
	step_size (array): maximum value to shift
	delta_t (float): delta time between tow points on the configuration (N/t_fin)
	noise (array): noise term of the equation
	metro_tol (float): maximum tolerance value to accept or reject a shift
"""
function acc_rej_move(config, N, chosen, step_size, delta_t, noise, metro_tol)
	start_resids = get_resids(config, delta_t, noise)
    shift_matrix = move_position(N, chosen, step_size)
	new_resids = get_resids(config + shift_matrix, delta_t, noise)
    exp_diff = exp.(new_resids - start_resids)
	
    checking = exp_diff[chosen-2] <= metro_tol # The acceptance rule must be evaluated in the position i-2 
												# because we are adjusting the second order derivative
    if checking
		return config + shift_matrix, 1
	else
		return config, 0
	end
end

"""
	main_here(mc_steps, step_size, N, t_fin, noise, error_tol, metro_tol, x0, v0)
Main loop for the Monte Carlo of the Next Guess Algorithm
ARGS
	
	mc_steps (int): Monte Carlo iterations
	step_size (array): maximum value to shift when moving
	time_told (int): lenght of the array containing the solution
	t_fin (int): lenght of the time interval to solve. Notice that delta_t = time_told/t_fin
	noise (array): noise term of the equation
	error_tol (float): maximum value accepted for the residuals of the whole system
	metro_tol (float): maximum tolerance value to accept or reject a shift
	x0 (float): initial position of the system
	v0 (float): initial velocity of the system
"""
function main_here(mc_steps, step_size, N, t_fin, noise, error_tol, metro_tol, x0, v0, sampling = false)
	delta_t = t_fin/N
	running_config = get_first_guess(N, delta_t, x0, v0)
	
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
			movement = acc_rej_move(running_config, N, k, step_size, delta_t, noise, metro_tol)
			running_config = movement[1]
			acc_rate += movement[2]
			trial_index += 1
		end
			
		# Update wrong points
		num_wrong = []
		current_res = get_resids(running_config, delta_t, noise)
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
			println("Solution Found in $i Steps")
			return running_config, time_config, time_resids
		end
		
		# interface data
		if i%(mc_steps*0.01) == 0
			println("Running:", " ", 100*i/mc_steps, "%, ", "Acceptance: ", acc_rate, "/", trial_index, ", Number Wrong: ", length(num_wrong))
		end
	end
	
	println("No Solution")
	return running_config, time_config, time_resids
end

######################### ######################### #########################
#########################       RUN ALGORITHM       #########################
######################### ######################### #########################
h = 0.5
final_time = 30
N = 500
noise_steps = 1

noise = get_noise(h, noise_steps*N, final_time, 1, DATA_PATH)

error_tol = 0.0001
mc_steps = 111000000
metro_tol = 1.000001
step_size = 0.001#0.008*final_time/time_steps
x0 = 0
v0 = 3
times = [i*final_time/N for i in 0:N-1]

num_soln = main_here(mc_steps, step_size, N, final_time, noise, error_tol, metro_tol, x0, v0)
write_data_hdf5(DATA_PATH * "2o_sde-h-0.5-noise$N.hdf5", (times, noise))
write_data_hdf5(DATA_PATH * "2o_sde-h-0.5-$N-v0$v0.hdf5", (times, num_soln[1]))