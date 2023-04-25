#= 
@authors: 17thsaint, davidsantiagoquevedo
Adapted from: https://github.com/17thSaint/finance-thesis
=#

module utils
export make_path, read_hdf5_data, write_data_hdf5

using HDF5

function make_path(h, which, dir = "")
    path = "$dir"*"fBM-h-$h-$which.hdf5"
    return path
end

function read_hdf5_data(h, which, dir = "", slice = false, len = 10000)
	path = make_path(h, which, dir)
	file = h5open(path, "r")
	data = [read(file["values"], "deets_t"), read(file["values"], "deets_v")]
	if slice
		data = [read(file["values"], "deets_t")[1:len], read(file["values"], "deets_v")[1:len]]
	end
	return data
end

function write_data_hdf5(path, data)
	println("Writing: $path")
	binary_file_pos = h5open(path, "w")
	create_group(binary_file_pos,"values")
	vals = binary_file_pos["values"]
    vals["deets_v"] = data[2]
	vals["deets_t"] = data[1]
	close(binary_file_pos)
	println("Data Added, file closed: $path")
end
end