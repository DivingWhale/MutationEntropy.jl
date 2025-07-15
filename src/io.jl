
function read_pae_for_round(data_path::String, mutation::AbstractString, round_val::Int)::Union{Matrix{Float64}, Nothing}
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

function read_paes(
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
                round_pae = read_pae_for_round(source_data_path, mutation, round_idx)
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

function get_variant_paes(
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
        round_pae = read_pae_for_round(source_data_path, mutation, round_idx)
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

function compute_single_distance_matrix(datadir::String, mutation::String, round_val::Int)
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

function get_dist_map(datadir::String, mutation::String, round_val::Int; cache_manager::Union{CacheManager, Nothing}=nothing)
    if cache_manager !== nothing
        cached_matrices = load_cached_data(cache_manager, mutation, :distance)
        if cached_matrices !== nothing && round_val <= length(cached_matrices)
            return cached_matrices[round_val]
        end
    end
    return compute_single_distance_matrix(datadir, mutation, round_val)
end
