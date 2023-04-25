#= 
@authors: davidsantiagoquevedo, 17thSaint
Adapted from: https://github.com/17thSaint/finance-thesis/blob/master/Codes/fracworking-lang.jl
=#

DATA_PATH = "data/"

include("../src/utils.jl")
include("../src/fBm_integration.jl")
using .utils, .fbm_integration, PyPlot
pygui(true)


function get_noise(h, n, t_fin, which = rand(1:10), dir ="", make_new = false)
	fBM = read_hdf5_data(h, which, dir, true, n+1)[2]
	if make_new
		fBM = frac_brown_wiki2(h, n, t_fin)[2]
	end
	return [t_fin*(fBM[i+1] - fBM[i])/n for i in 1:n]
end

"""
    get_first_guess(time_told)

First guess of the Monte Carlo Algorithm. 
For the Fractional Langevin equations, the first guess is the analytical solution of the integer order part of de differential equation.
However, for this first benchmark we use x = 0

ARGS
    time_told (int): 
"""
function get_first_guess(time_told)
	x0 = [0.0 for _ in 1:time_told]
    return x0
end

"""
    get_goft(h, config, delta_t, noise, noise_steps, lambda, gam)
Integer order part of the differential equation
Removed get_scale_inv_vals: Notice that in the original implementation they were returning 1.

"""
function get_goft(config, delta_t, noise, noise_steps)
	sliced_noise = [noise[i*noise_steps + 1] for i in 1:length(config) -2]
	velocity = [(config[i+1] - config[i])/delta_t for i in 1:length(config)-1]
	accel = [(velocity[i+1] - velocity[i])/delta_t for i in 1:length(velocity)-1]
	g_of_t = accel - sliced_noise
	return g_of_t, velocity, accel
end

function get_resids(config,delta_t,noise, noise_steps)
	g_stuff = get_goft(config, delta_t, noise, noise_steps)
	steps = length(config)
	right_hand_side = [0 for _ in 2:steps-1]
	return abs.(g_stuff[1] - right_hand_side)
end

function move_position(num_times, chosen, step_size)
	shift_matrix = [0.0 for i = 1:num_times]
	shift_matrix[chosen] += rand(-1:2:1)*rand(Float64)*step_size
	return shift_matrix
end

function acc_rej_move(config, num_times, chosen, step_size, delta_t, noise, noise_steps, top_val)
	start_resids = get_resids(config, delta_t, noise, noise_steps)
	
    shift_matrix = move_position(num_times, chosen, step_size)
	new_resids = get_resids(config + shift_matrix, delta_t, noise, noise_steps)
	
    exp_diff = exp.(new_resids - start_resids)
    checking = exp_diff[chosen-2] <= top_val
	
    if checking
		return config + shift_matrix, 1#, new_resids, start_resids, shift_matrix
	else
		return config, 0 #, new_resids, start_resids, shift_matrix
	end
end

function main_here(tol,steps,step_size, time_told, t_fin, noise, noise_steps, top_val, fixed)
	# getting first config from first guess
	running_config = append!([0.0,fixed],get_first_guess(time_told)[3:time_told])
	
    # only save configuration data for every 10 attempted movements
	samp_freq = 10

	time_config = fill(0.0,(time_told,Int(steps/samp_freq)))
	time_resids = fill(0.0,(time_told-2,Int(steps/samp_freq)))

	index = 1
	index2 = 0

	delta_t = t_fin/time_told
	acc_rate = 0

	num_wrong = [i+2 for i in 1:time_told-2]
	
	
	for i in 1:steps
		# each MC time step every time point has attempted move, except starting point
		#=
		for k in 3:time_told
			movement = acc_rej_move(running_config,h,time_told,k,step_size,delta_t,noise,noise_steps,top_val)
			running_config = movement[1]
			acc_rate += movement[2]
			#println(movement[2],", ",movement[3],movement[4],movement[5])
		end
		=#
		
		upper = 4
		if length(num_wrong) < 4
			upper = length(num_wrong)
		end
		
		for k in num_wrong[1:upper]#3 + num_correct:upper
			movement = acc_rej_move(running_config, time_told, k, step_size, delta_t, noise, 1, top_val)
			running_config = movement[1]
			acc_rate += movement[2]
			index2 += 1
		end
		num_wrong = []
		
		current_res = get_resids(running_config, delta_t, noise, noise_steps)
		
		# saving data every 10 steps
		if i%samp_freq == 0
			time_config[:,index] = [running_config[x] for x in 1:time_told]
			time_resids[:,index] = [current_res[y] for y in 1:time_told-2]
			index += 1
		end
		
		# if every time point has residuals less than tolerance then solution is found
		check_tol = [ current_res[j] < tol for j in 1:time_told-2 ]
		for i in 1:time_told-2
			if check_tol[i]
				
			else
				append!(num_wrong,[i+2])
			end
		end
		if all(check_tol)
			println("Solution Found in $i Steps")
			return running_config,time_config,time_resids,i
		end
		
		# interface data
		if i%(steps*0.01) == 0
			println("Running:"," ",100*i/steps,"%, ","Acceptance: ",acc_rate,"/",index2,", Number Wrong: ",length(num_wrong))
		end
	end
	
	println("No Solution")
	return time_config, time_resids
end


h = 0.5

final_time = 20
time_steps = 20
noise_steps = 1

noise = get_noise(h, noise_steps*time_steps, final_time, 1, DATA_PATH)
plot(noise)

times = [i*final_time/time_steps for i in 0:time_steps-1]


tol = 0.0001
mc_steps = 1000000
metro_val = 1.000001
step_size = 0.001#0.008*final_time/time_steps


num_soln = main_here(tol, mc_steps, step_size, time_steps, final_time, noise, noise_steps, metro_val, 0)

plot(num_soln[1])