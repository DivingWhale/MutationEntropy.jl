function read_pae_for_round_legacy(data_path::String, mutation::AbstractString, round_val::Int)::Union{Matrix{Float64}, Nothing}
    round_dir = joinpath(data_path, "$(mutation)_$(round_val)")
    !isdir(round_dir) && return nothing

    seed_dirs = sort(glob("seed-*_sample-0", round_dir))
    isempty(seed_dirs) && return nothing

    confidence_file = joinpath(first(seed_dirs), "confidences.json")
    !isfile(confidence_file) && return nothing

    try
        raw_data = JSON.parsefile(confidence_file)
        pae_data = raw_data["pae"]
        return reduce(hcat, [map(Float64, row) for row in pae_data])'
    catch e
        @warn "Failed to read or parse PAE data from $confidence_file: $e"
        return nothing
    end
end

function read_pae_from_path(pae_path::String)::Union{Matrix{Float64}, Nothing}
    !isdir(pae_path) && return nothing

    confidence_file = joinpath(pae_path, "confidences.json")
    !isfile(confidence_file) && return nothing

    try
        raw_data = JSON.parsefile(confidence_file)
        pae_data = raw_data["pae"]
        return reduce(hcat, [map(Float64, row) for row in pae_data])'
    catch e
        @warn "Failed to read or parse PAE data from $confidence_file: $e"
        return nothing
    end
end

function read_paes_legacy(
    task_file_path::String,
    source_data_path::String;
    verbose::Bool = false,
    num_rounds::Int = 20,
    cache_manager::Union{CacheManager,Nothing} = nothing,
)::Tuple{Dict{String,Vector{Matrix{Float64}}},CacheManager}
    if cache_manager === nothing
        cache_dir = joinpath(source_data_path, "cache")
        cache_manager = get_cache_manager(cache_dir)
    end

    verbose && println("Initializing multi-round PAE reading with cache.")

    mutations_list = read_mutations_from_file(task_file_path)
    paes = Dict{String,Vector{Matrix{Float64}}}()
    mutations_to_process = String[]

    for mutation in mutations_list
        cached_pae = load_cached_data(cache_manager, mutation, :pae)
        if cached_pae !== nothing
            paes[mutation] = cached_pae
            verbose && println("Loaded $mutation from cache.")
        else
            push!(mutations_to_process, mutation)
        end
    end

    verbose && println("Cache hits: $(length(paes))/$(length(mutations_list)). Processing $(length(mutations_to_process)) mutations.")

    if !isempty(mutations_to_process)
        results_channel = Channel{Tuple{String,Vector{Matrix{Float64}}}}(length(mutations_to_process))

        Threads.@threads for mutation in mutations_to_process
            pae_matrices = Vector{Matrix{Float64}}()
            for round_idx in 1:num_rounds
                round_pae = read_pae_for_round_legacy(source_data_path, mutation, round_idx)
                round_pae !== nothing && push!(pae_matrices, round_pae)
            end

            if !isempty(pae_matrices)
                put!(results_channel, (mutation, pae_matrices))
            elseif verbose
                @warn "No PAE data found for $mutation."
            end
        end
        close(results_channel)

        new_data = Dict{String,Vector{Matrix{Float64}}}()
        for (mutation, matrices) in results_channel
            paes[mutation] = matrices
            new_data[mutation] = matrices
        end

        if !isempty(new_data)
            verbose && println("Saving $(length(new_data)) new items to cache.")
            for (mutation, matrices) in new_data
                save_cached_data(cache_manager, mutation, :pae, matrices)
            end
        end
    end

    verbose && println("Finished processing. Total mutations with PAE data: $(length(paes)).")
    return paes, cache_manager
end

function get_variant_paes_legacy(
    mutation::String,
    source_data_path::String;
    verbose::Bool = false,
    num_rounds::Int = 20,
    cache_manager::Union{CacheManager,Nothing} = nothing,
)::Tuple{Union{Vector{Matrix{Float64}}, Nothing}, CacheManager}
    if cache_manager === nothing
        cache_dir = joinpath(source_data_path, "cache")
        cache_manager = get_cache_manager(cache_dir)
    end

    cached_pae = load_cached_data(cache_manager, mutation, :pae)
    if cached_pae !== nothing
        verbose && println("Loaded $mutation from cache.")
        return cached_pae, cache_manager
    end

    verbose && println("Cache miss for $mutation. Reading from source.")
    pae_matrices = Vector{Matrix{Float64}}()
    for round_idx in 1:num_rounds
        round_pae = read_pae_for_round_legacy(source_data_path, mutation, round_idx)
        round_pae !== nothing && push!(pae_matrices, round_pae)
    end

    if isempty(pae_matrices)
        verbose && @warn "No PAE data found for $mutation."
        return nothing, cache_manager
    end

    verbose && println("Saving $mutation to cache.")
    save_cached_data(cache_manager, mutation, :pae, pae_matrices)
    
    return pae_matrices, cache_manager
end

function get_variant_paes(
    mutation::String,
    source_data_path::String;
    verbose::Bool = false,
    num_rounds::Int = 20, # num_rounds can be used as an upper limit for seeds
    cache_manager::Union{CacheManager,Nothing} = nothing,
)::Tuple{Union{Vector{Matrix{Float64}}, Nothing}, CacheManager}
    if cache_manager === nothing
        cache_dir = joinpath(source_data_path, "cache")
        cache_manager = get_cache_manager(cache_dir)
    end

    cached_pae = load_cached_data(cache_manager, mutation, :pae)
    if cached_pae !== nothing
        verbose && println("Loaded $mutation from cache.")
        return cached_pae, cache_manager
    end

    verbose && println("Cache miss for $mutation. Reading from source.")
    mutation_dir = joinpath(source_data_path, mutation)
    if !isdir(mutation_dir)
        verbose && @warn "Mutation directory not found: $mutation_dir"
        return nothing, cache_manager
    end

    best_samples = get_best_samples(mutation_dir)
    if isempty(best_samples)
        verbose && @warn "No best samples found from ranking_scores.csv for $mutation."
        return nothing, cache_manager
    end

    pae_matrices = Vector{Matrix{Float64}}()
    # Sort seeds to maintain order
    sorted_seeds = sort(collect(keys(best_samples)))

    for seed_idx in sorted_seeds
        # Respect num_rounds as an upper limit
        seed_idx > num_rounds && break
        
        sample_idx = best_samples[seed_idx]
        pae_path = joinpath(mutation_dir, "seed-$(seed_idx)_sample-$(sample_idx)")
        
        round_pae = read_pae_from_path(pae_path)
        round_pae !== nothing && push!(pae_matrices, round_pae)
    end

    if isempty(pae_matrices)
        verbose && @warn "No PAE data successfully read for $mutation."
        return nothing, cache_manager
    end

    verbose && println("Saving $mutation to cache.")
    save_cached_data(cache_manager, mutation, :pae, pae_matrices)
    
    return pae_matrices, cache_manager
end

function compute_single_distance_matrix_legacy(datadir::String, mutation::String, round_val::Int)
    base_dir = joinpath(datadir, "$(mutation)_$(round_val)")
    subfolder_pattern = "seed-*_sample-0"
    
    if !isdir(base_dir)
        error("Base directory not found: $base_dir")
    end
    
    subfolders = glob(subfolder_pattern, base_dir)
    subfolders = [joinpath(base_dir, sf) for sf in subfolders]
    
    if isempty(subfolders)
        error("No matching folder found: $subfolder_pattern in $base_dir")
    end
    
    subfolder = subfolders[1]
    model_pdb = joinpath(subfolder, "model.pdb")
    
    _ensure_pdb_from_cif_if_needed(subfolder, model_pdb)

    d_matrix_file = joinpath(subfolder, "d_matrix.jld2")
    
    local d_matrix
    if !isfile(d_matrix_file)
        structure = read(model_pdb, PDBFormat)
        d_matrix = ca_distMap(structure).data
        
        try
            jldsave(d_matrix_file; d_matrix=d_matrix)
        catch e
            @warn "Failed to save distance matrix to $d_matrix_file: $e"
        end
    else
        try
            d_matrix = load(d_matrix_file, "d_matrix")
        catch e
            @warn "Failed to load cached distance matrix from $d_matrix_file: $e. Recomputing..."
            structure = read(model_pdb, PDBFormat)
            d_matrix = ca_distMap(structure).data
            try
                jldsave(d_matrix_file; d_matrix=d_matrix)
            catch save_err
                @warn "Failed to save recomputed distance matrix to $d_matrix_file: $save_err"
            end
        end
    end
    
    return d_matrix
end

function get_dist_map_legacy(datadir::String, mutation::String, round_val::Int; cache_manager::Union{CacheManager, Nothing}=nothing)
    if cache_manager !== nothing
        cached_matrices = load_cached_data(cache_manager, mutation, :distance)
        if cached_matrices !== nothing && round_val <= length(cached_matrices)
            return cached_matrices[round_val]
        end
    end
    return compute_single_distance_matrix_legacy(datadir, mutation, round_val)
end

function get_dist_maps_legacy(
    mutation::String,
    source_data_path::String;
    verbose::Bool = false,
    num_rounds::Int = 20,
    cache_manager::Union{CacheManager,Nothing} = nothing,
)::Tuple{Union{Vector{Matrix{Float64}}, Nothing}, CacheManager}
    if cache_manager === nothing
        cache_dir = joinpath(source_data_path, "cache")
        cache_manager = get_cache_manager(cache_dir)
    end

    cached_dist = load_cached_data(cache_manager, mutation, :distance)
    if cached_dist !== nothing
        verbose && println("Loaded distance matrices for $mutation from cache.")
        return cached_dist, cache_manager
    end

    verbose && println("Cache miss for distance matrices of $mutation. Reading from source.")
    dist_matrices = Vector{Matrix{Float64}}()
    for round_idx in 1:num_rounds
        dist_matrix = compute_single_distance_matrix_legacy(source_data_path, mutation, round_idx)
        if dist_matrix !== nothing
            push!(dist_matrices, dist_matrix)
        end
    end

    if isempty(dist_matrices)
        verbose && @warn "No distance matrices found for $mutation."
        return nothing, cache_manager
    end

    verbose && println("Saving distance matrices for $mutation to cache.")
    save_cached_data(cache_manager, mutation, :distance, dist_matrices)
    
    return dist_matrices, cache_manager
end

function compute_dist_matrix_from_path(model_dir::String)
    if !isdir(model_dir)
        @warn "Model directory not found: $model_dir"
        return nothing
    end
    
    model_pdb = joinpath(model_dir, "model.pdb")
    d_matrix_file = joinpath(model_dir, "d_matrix.jld2")

    # First, check for a pre-computed distance matrix
    if isfile(d_matrix_file)
        try
            return load(d_matrix_file, "d_matrix")
        catch e
            @warn "Failed to load cached distance matrix from $d_matrix_file: $e. Recomputing..."
        end
    end

    # If not found or failed to load, compute it from the PDB file
    _ensure_pdb_from_cif_if_needed(model_dir, model_pdb)
    if !isfile(model_pdb)
        @warn "PDB file not found after ensuring conversion: $model_pdb"
        return nothing
    end

    try
        structure = read(model_pdb, PDBFormat)
        d_matrix = ca_distMap(structure).data
        
        # Save the newly computed matrix
        try
            jldsave(d_matrix_file; d_matrix=d_matrix)
        catch e
            @warn "Failed to save distance matrix to $d_matrix_file: $e"
        end
        
        return d_matrix
    catch e
        @error "Failed to read PDB or compute distance matrix for $model_pdb: $e"
        return nothing
    end
end

function get_dist_maps(
    mutation::String,
    source_data_path::String;
    verbose::Bool = false,
    num_rounds::Int = 20, # Acts as an upper limit for seeds
    cache_manager::Union{CacheManager,Nothing} = nothing,
)::Tuple{Union{Vector{Matrix{Float64}}, Nothing}, CacheManager}
    if cache_manager === nothing
        cache_dir = joinpath(source_data_path, "cache")
        cache_manager = get_cache_manager(cache_dir)
    end

    # Check for cached data for the entire mutation
    cached_dist_maps = load_cached_data(cache_manager, mutation, :distance)
    if cached_dist_maps !== nothing
        verbose && println("Loaded distance maps for $mutation from cache.")
        return cached_dist_maps, cache_manager
    end

    verbose && println("Cache miss for distance maps of $mutation. Reading from source.")
    mutation_dir = joinpath(source_data_path, mutation)
    if !isdir(mutation_dir)
        verbose && @warn "Mutation directory not found: $mutation_dir"
        return nothing, cache_manager
    end

    best_samples = get_best_samples(mutation_dir)
    if isempty(best_samples)
        verbose && @warn "No best samples found from ranking_scores.csv for $mutation."
        return nothing, cache_manager
    end

    dist_matrices = Vector{Matrix{Float64}}()
    sorted_seeds = sort(collect(keys(best_samples)))

    for seed_idx in sorted_seeds
        seed_idx > num_rounds && break
        
        sample_idx = best_samples[seed_idx]
        model_path = joinpath(mutation_dir, "seed-$(seed_idx)_sample-$(sample_idx)")
        
        dist_matrix = compute_dist_matrix_from_path(model_path)
        dist_matrix !== nothing && push!(dist_matrices, dist_matrix)
    end

    if isempty(dist_matrices)
        verbose && @warn "No distance matrices successfully computed for $mutation."
        return nothing, cache_manager
    end

    verbose && println("Saving distance maps for $mutation to cache.")
    save_cached_data(cache_manager, mutation, :distance, dist_matrices)
    
    return dist_matrices, cache_manager
end


function get_best_samples(mutation_dir::String)::Dict{Int, Int}
    ranking_file = joinpath(mutation_dir, "ranking_scores.csv")
    !isfile(ranking_file) && return Dict{Int, Int}()

    try
        df = CSV.read(ranking_file, DataFrame)
        best_samples = Dict{Int, Int}()
        
        # Group by seed and find the sample with the highest ranking_score for each seed
        gdf = DataFrames.groupby(df, :seed)
        for sub_df in gdf
            if !isempty(sub_df)
                best_row = sub_df[argmax(sub_df.ranking_score), :]
                best_samples[best_row.seed] = best_row.sample
            end
        end
        return best_samples
    catch e
        @warn "Failed to read or process ranking scores from $ranking_file: $e"
        return Dict{Int, Int}()
    end
end
