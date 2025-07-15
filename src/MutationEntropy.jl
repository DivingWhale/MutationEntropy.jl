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
using SHA
using Printf

# Core functionality exports
export read_coordinates, read_xvg, calculate_ddgs, ΔΔS
export get_residue_calpha_b_factor

# New structured data exports
export MutationData, EntropyParams

# New cache system exports
export CacheManager, get_cache_manager
export has_cached_data, load_cached_data, save_cached_data
export preload_cached_data, clear_cache, clear_all_cache
export get_cache_stats, print_cache_stats, get_cache_size_info

include("types.jl")
include("utils.jl")
include("cache_manager.jl")
include("io.jl")
include("entropy.jl")
include("ddg.jl")
include("mutation_effect.jl")
include("plot_functions.jl")

end
