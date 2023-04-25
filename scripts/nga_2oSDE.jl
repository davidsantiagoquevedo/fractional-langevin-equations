#= 
@authors: davidsantiagoquevedo, 17thSaint
Adapted from: https://github.com/17thSaint/finance-thesis
=#

DATA_PATH = "data/"

include("../src/utils.jl")
include("../src/fBm_integration.jl")
using .utils, .fbm_integration

