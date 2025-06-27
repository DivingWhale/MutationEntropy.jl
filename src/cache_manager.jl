"""
Cache management system for MutationEntropy.jl

This module provides efficient caching for PAE and distance matrices.
All data is stored per-mutation with all rounds in a single file.
No backward compatibility or migration code is included.
"""

using JLD2
using Printf
import Base: isdir, mkpath

# Cache directory structure
const DEFAULT_CACHE_DIR = "cache"
const PAE_CACHE_SUBDIR = "pae"
const DIST_CACHE_SUBDIR = "distance"

"""
CacheManager handles caching of PAE and distance matrices.
All matrices are stored per-mutation with all rounds in a single file.
"""
mutable struct CacheManager
    cache_dir::String
    pae_cache_dir::String
    dist_cache_dir::String
    stats::Dict{String, Int}
    
    function CacheManager(cache_dir::String = DEFAULT_CACHE_DIR)
        pae_dir = joinpath(cache_dir, PAE_CACHE_SUBDIR)
        dist_dir = joinpath(cache_dir, DIST_CACHE_SUBDIR)
        
        # Create directories if they don't exist
        for dir in [cache_dir, pae_dir, dist_dir]
            if !isdir(dir)
                mkpath(dir)
            end
        end
        
        stats = Dict("pae_hits" => 0, "pae_misses" => 0, 
                    "dist_hits" => 0, "dist_misses" => 0,
                    "pae_writes" => 0, "dist_writes" => 0)
        
        new(cache_dir, pae_dir, dist_dir, stats)
    end
end

"""
Get cache file path for a given mutation and data type.
All files use the '_all_rounds' suffix to indicate they contain all rounds.
"""
function get_cache_path(cache_manager::CacheManager, mutation::String, data_type::Symbol)
    cache_key = "$(mutation)_all_rounds"
    
    if data_type == :pae
        return joinpath(cache_manager.pae_cache_dir, "$(cache_key).jld2")
    elseif data_type == :distance
        return joinpath(cache_manager.dist_cache_dir, "$(cache_key).jld2")
    else
        error("Unknown data type: $data_type. Must be :pae or :distance")
    end
end

"""
Check if cached data exists for a mutation.
"""
function has_cached_data(cache_manager::CacheManager, mutation::String, data_type::Symbol)
    cache_path = get_cache_path(cache_manager, mutation, data_type)
    return isfile(cache_path)
end

"""
Load cached data for a mutation (all rounds).
Returns a vector of matrices, one per round.
"""
function load_cached_data(cache_manager::CacheManager, mutation::String, data_type::Symbol)
    cache_path = get_cache_path(cache_manager, mutation, data_type)
    
    if !isfile(cache_path)
        if data_type == :pae
            cache_manager.stats["pae_misses"] += 1
        else
            cache_manager.stats["dist_misses"] += 1
        end
        return nothing
    end
    
    try
        data = load_object(cache_path)
        if data_type == :pae
            cache_manager.stats["pae_hits"] += 1
        else
            cache_manager.stats["dist_hits"] += 1
        end
        return data
    catch e
        @warn "Failed to load cache file $cache_path: $e"
        return nothing
    end
end

"""
Save data to cache for a mutation (all rounds).
Data should be a vector of matrices, one per round.
"""
function save_cached_data(cache_manager::CacheManager, mutation::String, data_type::Symbol, data::Vector)
    cache_path = get_cache_path(cache_manager, mutation, data_type)
    
    try
        # Create directory if it doesn't exist
        cache_dir = dirname(cache_path)
        if !isdir(cache_dir)
            mkpath(cache_dir)
        end
        
        save_object(cache_path, data)
        
        if data_type == :pae
            cache_manager.stats["pae_writes"] += 1
        else
            cache_manager.stats["dist_writes"] += 1
        end
        
        @info "Cached $data_type data for $mutation ($(length(data)) rounds)"
        return true
    catch e
        @warn "Failed to save cache file $cache_path: $e"
        return false
    end
end

"""
Preload multiple mutations into memory for batch processing.
Returns a dictionary mapping mutation -> data.
"""
function preload_cached_data(cache_manager::CacheManager, mutations::Vector{String}, data_type::Symbol)
    preloaded = Dict{String, Vector}()
    
    for mutation in mutations
        data = load_cached_data(cache_manager, mutation, data_type)
        if data !== nothing
            preloaded[mutation] = data
        end
    end
    
    @info "Preloaded $data_type data for $(length(preloaded))/$(length(mutations)) mutations"
    return preloaded
end

"""
Clear all cache data for a specific data type.
"""
function clear_cache(cache_manager::CacheManager, data_type::Symbol)
    cache_dir = data_type == :pae ? cache_manager.pae_cache_dir : cache_manager.dist_cache_dir
    
    if isdir(cache_dir)
        try
            rm(cache_dir, recursive=true)
            mkpath(cache_dir)
            @info "Cleared $data_type cache directory"
        catch e
            @warn "Failed to clear $data_type cache: $e"
        end
    end
end

"""
Clear all cache data.
"""
function clear_all_cache(cache_manager::CacheManager)
    clear_cache(cache_manager, :pae)
    clear_cache(cache_manager, :distance)
    
    # Reset stats
    for key in keys(cache_manager.stats)
        cache_manager.stats[key] = 0
    end
end

"""
Get cache statistics.
"""
function get_cache_stats(cache_manager::CacheManager)
    stats = Dict{String, Any}()
    
    # Copy integer stats
    for (key, value) in cache_manager.stats
        stats[key] = value
    end
    
    # Calculate hit rates
    pae_total = stats["pae_hits"] + stats["pae_misses"]
    dist_total = stats["dist_hits"] + stats["dist_misses"]
    
    if pae_total > 0
        stats["pae_hit_rate"] = stats["pae_hits"] / pae_total
    else
        stats["pae_hit_rate"] = 0.0
    end
    
    if dist_total > 0
        stats["dist_hit_rate"] = stats["dist_hits"] / dist_total
    else
        stats["dist_hit_rate"] = 0.0
    end
    
    return stats
end

"""
Print cache statistics in a readable format.
"""
function print_cache_stats(cache_manager::CacheManager)
    stats = get_cache_stats(cache_manager)
    
    println("=== Cache Statistics ===")
    println("PAE Cache:")
    println("  Hits: $(stats["pae_hits"])")
    println("  Misses: $(stats["pae_misses"])")
    println("  Hit Rate: $(@sprintf("%.1f", stats["pae_hit_rate"] * 100))%")
    println("  Writes: $(stats["pae_writes"])")
    
    println("Distance Cache:")
    println("  Hits: $(stats["dist_hits"])")
    println("  Misses: $(stats["dist_misses"])")
    println("  Hit Rate: $(@sprintf("%.1f", stats["dist_hit_rate"] * 100))%")
    println("  Writes: $(stats["dist_writes"])")
    println("========================")
end

"""
Get cache size information.
"""
function get_cache_size_info(cache_manager::CacheManager)
    pae_files = 0
    dist_files = 0
    pae_size = 0
    dist_size = 0
    
    if isdir(cache_manager.pae_cache_dir)
        for file in readdir(cache_manager.pae_cache_dir)
            if endswith(file, ".jld2")
                pae_files += 1
                path = joinpath(cache_manager.pae_cache_dir, file)
                pae_size += stat(path).size
            end
        end
    end
    
    if isdir(cache_manager.dist_cache_dir)
        for file in readdir(cache_manager.dist_cache_dir)
            if endswith(file, ".jld2")
                dist_files += 1
                path = joinpath(cache_manager.dist_cache_dir, file)
                dist_size += stat(path).size
            end
        end
    end
    
    return Dict(
        "pae_files" => pae_files,
        "dist_files" => dist_files,
        "pae_size_mb" => pae_size / (1024^2),
        "dist_size_mb" => dist_size / (1024^2),
        "total_size_mb" => (pae_size + dist_size) / (1024^2)
    )
end

# Global cache manager instance
const GLOBAL_CACHE_MANAGER = Ref{Union{CacheManager, Nothing}}(nothing)

"""
Get or create the global cache manager instance.
"""
function get_cache_manager(cache_dir::String = DEFAULT_CACHE_DIR)
    if GLOBAL_CACHE_MANAGER[] === nothing || GLOBAL_CACHE_MANAGER[].cache_dir != cache_dir
        GLOBAL_CACHE_MANAGER[] = CacheManager(cache_dir)
    end
    return GLOBAL_CACHE_MANAGER[]
end

# Export public interface
export CacheManager, get_cache_manager
export has_cached_data, load_cached_data, save_cached_data
export preload_cached_data, clear_cache, clear_all_cache
export get_cache_stats, print_cache_stats, get_cache_size_info
