module MutationEntropy

using BioStructures, LinearAlgebra
using CairoMakie
using Glob
using JSON
using IterTools
using Statistics
using DataFrames

export read_coordinates, read_xvg, plot_msf, calculate_ddgs

include("entropy.jl")
include("visualization.jl")
include("ddg.jl")

end
