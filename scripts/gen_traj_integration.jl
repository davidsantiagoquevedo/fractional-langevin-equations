#=
@author: davidsantiagoquevedo
Adapted from: https://github.com/17thSaint/finance-thesis
=#

DATA_PATH = "data/fbm/"

include("../src/utils.jl")
include("../src/fBm_integration.jl")
using .utils, .fbm_integration

h = parse(Float64,ARGS[1])
n = parse(Int64,ARGS[2])
T = parse(Int64,ARGS[3])
trajectory_i = parse(Int64,ARGS[4])
trajectory_f = parse(Int64,ARGS[5])

for j in trajectory_i:trajectory_f
    path = "$DATA_PATH"*"fBM-h-$h-$j.hdf5"
	trajectory = frac_brown_wiki2(h,n,T)
	write_data_hdf5(path,trajectory)
end