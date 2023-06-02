#= 
@authors: 17thsaint, davidsantiagoquevedo
Adapted from: https://github.com/17thSaint/finance-thesis
=#

module utils
export read_hdf5, write_data_hdf5

using HDF5

function read_hdf5(path)
	file = h5open(path, "r")
	data = [read(file["values"], "deets_t"), read(file["values"], "deets_v")]
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