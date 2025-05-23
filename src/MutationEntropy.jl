module MutationEntropy

using BioStructures, LinearAlgebra
using CairoMakie
using Glob
using JSON
using IterTools
using Statistics
using DataFrames
using JLD2
using Distributed

export read_coordinates, read_xvg, plot_msf, calculate_ddgs

include("entropy.jl")
include("visualization.jl")
include("ddg.jl")
include("mutation_effect.jl")
include("plot_functions.jl")

end
