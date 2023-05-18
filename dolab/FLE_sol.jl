using HDF5, Cubature, MittagLeffler, SpecialFunctions, Plots

function read_hdf5_data(h,count,slice=false,len=10000)
	file = h5open("data/fBM-h-$h-$count.hdf5","r")
	data = [read(file["values"],"deets_t"), read(file["values"],"deets_v")]
	if slice
		data = [read(file["values"],"deets_t")[1:len], read(file["values"],"deets_v")[1:len]]
	end
	return data
end


# gets noise from previous stored data file or makes new if asked
function noise(h,n,t_fin,which=rand(1:20),make_new=false)
	fBM = read_hdf5_data(h,which,true,n+1)[2]
	if make_new
		fBM = frac_brown_wiki2(h,n,t_fin)[2]
	end
	return [t_fin*(fBM[i+1]-fBM[i])/n for i in 1:n]
end

function eta(gam,h,kBT)
	#kBT = 1.0
	return sqrt(2*gam*kBT*gamma(1.5-h)*gamma(0.5+h)/(gamma(2*h)*gamma(2-2*h)))
end

function lang_soln(h,t_steps,noise_steps,gam,m,t_fin,v0,which=rand(1:20))
	times = [i*t_fin/t_steps for i in 0:t_steps]
	term_one = [0.0 for i in 1:t_steps]
	term_two = [0.0 for i in 1:t_steps]
	c_eta = eta(gam,h,1.0)
	noise_term = c_eta.*noise(h,Int(t_steps*noise_steps),t_fin,which)
	for i in 1:t_steps
		if i%(0.05*t_steps) == 0
			println("H=",h,", ",100*i/t_steps,"%",", Lang")
		end
		noise_times = append!([0.0],[ l*t_fin/t_steps + j*t_fin/(t_steps*noise_steps) for l in 0:i-1 for j in 1:noise_steps])
		for k in 1:length(noise_times)-1
			left_time = noise_times[k]
			right_time = noise_times[k+1]
			
			mitlef_left = mittleff(2*h,2,(-gam/m)*(left_time^(2*h)))
			mitlef_right = mittleff(2*h,2,(-gam/m)*(right_time^(2*h)))
				
			if v0 != 0.0 && k == length(noise_times)-1
				term_two[i] += v0*right_time*mitlef_right
			end
			
			if k == length(noise_times)-1
				noise_left = noise_term[1]
			else
				noise_left = noise_term[length(noise_times)-1-k]
			end
			noise_right = noise_term[length(noise_times)-k]
				
			term_one[i] += 0.5*(right_time-left_time)*c_eta*(noise_left*left_time*mitlef_left + noise_right*right_time*mitlef_right )
		end
		
	end
	full_position = append!([0.0],term_one./m + term_two)
	return times,full_position
end



time_steps = 100
gam = 1.0
mass = 1.0
final_time = 10
v0 = 0.0
solns = [lang_soln(i,time_steps,10,gam,mass,final_time,v0,2) for i in [0.3]]
v0 = 0.0002
solns1 = [lang_soln(i,time_steps,10,gam,mass,final_time,v0,2) for i in [0.3]]

plot(solns, marker =:circle)
plot!(solns1)