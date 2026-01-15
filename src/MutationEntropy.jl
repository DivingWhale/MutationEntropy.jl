module MutationEntropy

using BioStructures, LinearAlgebra
using CairoMakie
using SwarmMakie
using Glob
using JSON
using IterTools
using Statistics
using DataFrames
using JLD2
using Distributed
using SHA
using Printf
using GLM
using CSV
using TOML
using Dates
using Colors

# Core functionality exports
export read_coordinates, read_xvg, calculate_ddgs, ΔΔS, read_mutations_from_file, parse_mutation_position, get_variant_paes, get_dist_map, find_residues_within_distance, read_ddg_rosetta, read_ddg_foldx, read_ddg_pythia, read_ddg_stabilityoracle, read_experimental_csv, read_pythia_txt
export perform_correlation_analysis, process_entropy_data, plot_correction_scatter
export get_residue_calpha_b_factor, write_b_factors_to_pdb
export get_variant_paes_legacy, get_dist_maps_legacy  # Legacy data structure support
export calculate_mse_per_position, plot_swarm_by_position

# PCC analysis exports
export analyze_pcc_files, analyze_fixed_rho, load_pcc_results, read_correlation_data
export find_best_pcc_parameters, save_pcc_results

# PCC plotting exports
export plot_pcc_vs_A, get_baseline_correlation

# Configuration exports
export ProteinConfig, PredictorConfig, load_protein_config, get_predictor_config
export create_protein_config, config_from_legacy, generate_config_toml
export DEFAULT_ANALYSIS_PARAMS


# New structured data exports
export MutationData, EntropyParams

# New cache system exports
export CacheManager, get_cache_manager
export has_cached_data, load_cached_data, save_cached_data
export preload_cached_data, clear_cache, clear_all_cache
export get_cache_stats, print_cache_stats, get_cache_size_info

# Multi-mutation support exports
export is_multi_mutation, parse_multi_mutation, parse_mutation_positions

include("types.jl")
include("utils.jl")
include("cache_manager.jl")
include("io.jl")
include("entropy.jl")
include("ddg.jl")
include("mutation_effect.jl")
include("plot_functions.jl")
include("pcc_analysis.jl")
include("pcc_plotting.jl")
include("config.jl")

include("analysis.jl")

end
