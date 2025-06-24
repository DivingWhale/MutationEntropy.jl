Φ(r, η, κ) = exp(-(r/η)^κ)

function Γ(coordinates::AbstractVector{Vector{Float64}}, η::Int, κ::Int)
    natoms = length(coordinates)
    Γ = zeros(natoms, natoms)
    for i = 1:natoms, j = i+1:natoms
        Γ[i, j] = -Φ(norm(coordinates[i] - coordinates[j]), η, κ)
    end
    Γ += Γ'
    for i = 1:natoms
        Γ[i, i] = -sum(view(Γ, :, i))
    end
    return Γ
end

function compute_Γ(coordinates::AbstractVector{Vector{Float64}})
    Γ1 = Γ(coordinates, 20, 7)
    Γ2 = Γ(coordinates, 13, 10)
    return Γ1 + Γ2
end

function read_coordinates(filename::String)
    WT_structure = read(filename, PDBFormat)
    CA_atoms = collectatoms(WT_structure, sel"name CA")
    coordinates = getfield.(CA_atoms, :coords)
    return coordinates
end

function read_xvg(filename::String)
    result = Tuple{Int, Float64}[]
    for line in eachline(filename)
        if startswith(line, "#") || startswith(line, "@") || isempty(line)
            continue
        end
        parts = split(line)
        push!(result, (parse(Int, parts[1]), parse(Float64, parts[2])))
    end
    return result
end

function msf(Γ::AbstractMatrix)
    return diag(pinv(Γ))
end

function read_single_pae(data_path::String, mutation::AbstractString, round_val::Int)
    round_dir = joinpath(data_path, "$(mutation)_$(round_val)")
    
    if !isdir(round_dir)
        return nothing
    end
    
    # Find seed directories in this round and sort to get consistent ordering
    seed_dirs = sort(glob("seed-*_sample-0", round_dir))
    
    if isempty(seed_dirs)
        return nothing
    end
    
    # Read only the first seed (top1)
    first_seed_dir = seed_dirs[1]
    full_seed_path = joinpath(round_dir, first_seed_dir)
    confidence_file = joinpath(full_seed_path, "confidences.json")
    
    if !isfile(confidence_file)
        return nothing
    end
    
    try
        raw_data = JSON.parsefile(confidence_file)
        pae_any = raw_data["pae"]
        pae = map(x -> map(Float64, x), pae_any)
        pae = hcat(pae...)'
        
        return pae
    catch e
        @warn "Failed to read PAE data from $confidence_file: $e"
        return nothing
    end
end

function read_paes(task_file_path::String, source_data_path::String; verbose::Bool=false, num_rounds::Int=20)
    CACHE_DIR = joinpath(source_data_path, ".pae_cache")
    if verbose
        println("Initializing multi-round PAE reading. Cache directory: $(abspath(CACHE_DIR))")
    end

    try
        mkpath(CACHE_DIR)
        if verbose
            println("Cache directory ensured at $(abspath(CACHE_DIR))")
        end
    catch e
        @warn "Could not create cache directory $CACHE_DIR: $e. Caching may not work as expected."
    end

    mutations_list = read_mutations_from_file(task_file_path)
    num_mutations = length(mutations_list)
    if verbose
        println("Found $(num_mutations) mutations to process from $task_file_path.")
    end
    
    results_channel = Channel{Tuple{String, Union{Vector{Matrix{Float64}}, Nothing}}}(num_mutations)

    if verbose
        println("Starting parallel processing of mutations...")
    end
    processed_count_atomic = Threads.Atomic{Int}(0) # Atomic counter for progress

    Threads.@threads for mutation_original_case in mutations_list
        mutation = string(strip(mutation_original_case))
        
        current_processed_count = Threads.atomic_add!(processed_count_atomic, 1)
        if verbose
            println("Processing $(current_processed_count)/$(num_mutations): $mutation")
        end

        cache_file_name = "$(mutation)_multiround.jld2"
        cache_file_path = joinpath(CACHE_DIR, cache_file_name)
        
        loaded_pae = nothing

        if isfile(cache_file_path)
            try
                loaded_pae = JLD2.load_object(cache_file_path)
            catch err
                @warn "Failed to load $mutation from cache file $cache_file_path: $err. Recomputing."
            end
        end

        if loaded_pae === nothing # Cache miss or failed to load
            # Read PAE data for multiple rounds
            pae_matrices = Vector{Matrix{Float64}}()
            
            for round_idx in 1:num_rounds
                round_pae = read_single_pae(source_data_path, mutation, round_idx)
                if round_pae !== nothing
                    push!(pae_matrices, round_pae)
                end
            end
            
            loaded_pae = pae_matrices
            
            if !isempty(loaded_pae)
                try
                    JLD2.save_object(cache_file_path, loaded_pae)
                catch err
                    @warn "Failed to save $mutation to cache file $cache_file_path: $err."
                end
            else
                if verbose
                    println("No PAE data found for $mutation from source $source_data_path.") # Kept as it's an important specific outcome
                end
                loaded_pae = nothing
            end
        end
        
        put!(results_channel, (mutation, loaded_pae))
    end
    
    close(results_channel)
    if verbose
        println("Finished parallel processing. Aggregating results...")
    end

    paes = Dict{String, Vector{Matrix{Float64}}}()
    aggregated_count = 0
    for (mut, pae_matrices) in results_channel
        if pae_matrices !== nothing && !isempty(pae_matrices)
            paes[mut] = pae_matrices
        else
            if verbose
                println("Warning: No PAE data for mutation $mut after processing.")
            end
        end
        aggregated_count += 1
    end
    
    if verbose
        println("Finished aggregation. Processed $(num_mutations) attempted mutations. $(length(paes)) mutations have PAE data.")
    end
    
    return paes
end

# Configuration struct for entropy calculations
"""
    EntropyConfig

Configuration structure for entropy calculations.

# Fields
- `datadir::String`: Directory containing structure data
- `wt_identifier::String`: Identifier for wild-type structure (default: "WT")
- `wt_dist_matrices::Union{Vector{Matrix{Float64}}, Nothing}`: Pre-computed WT distance matrices for all rounds (optional)
- `offset::Int64`: Offset between protein position numbers and matrix indices (default: 0)
"""
struct EntropyConfig
    datadir::String
    wt_identifier::String
    wt_dist_matrices::Union{Vector{Matrix{Float64}}, Nothing}
    offset::Int64
    
    # Constructor with defaults
    function EntropyConfig(datadir::String; 
                          wt_identifier::String="WT", 
                          wt_dist_matrices::Union{Vector{Matrix{Float64}}, Nothing}=nothing,
                          offset::Int64=0)
        new(datadir, wt_identifier, wt_dist_matrices, offset)
    end
end

function ΔΔS(position::Int64, rho::Float64, α::Float64, Γ::Matrix{Float64}, 
             PAE_mut_rounds::Vector{Matrix{Float64}}, PAE_wt_rounds::Vector{Matrix{Float64}}, 
             mutation::AbstractString, config::EntropyConfig)
    matrix_idx = position - config.offset
    indices = findall(x -> abs(x) > 1e-5, Γ[matrix_idx, :])

    # Initialize arrays to store values for each round
    num_rounds = min(length(PAE_mut_rounds), length(PAE_wt_rounds))
    if num_rounds == 0
        @warn "No matching rounds found for mutation $mutation"
        return 0.0
    end

    # Calculate terms for each round and collect them
    mut_terms = zeros(Float64, num_rounds, length(indices))
    wt_terms = zeros(Float64, num_rounds, length(indices))
    
    for round_idx in 1:num_rounds
        PAE_mut = PAE_mut_rounds[round_idx]
        PAE_wt = PAE_wt_rounds[round_idx]
        
        # Get distance matrices for this round
        local dist_mut, dist_wt
        
        # For mutant, read distance matrix for this round
        try
            dist_mut = get_dist_map(config.datadir, String(mutation), round_idx)
        catch e
            @warn "Could not read mutant distance matrix for $mutation round $round_idx: $e. Using identity matrix."
            dist_mut = Matrix{Float64}(I, size(PAE_mut))
        end
        
        # For WT, use pre-computed matrices if available, otherwise read from disk
        if config.wt_dist_matrices !== nothing && round_idx <= length(config.wt_dist_matrices)
            # Use pre-computed matrix for this round
            dist_wt = config.wt_dist_matrices[round_idx]
        else
            try
                dist_wt = get_dist_map(config.datadir, config.wt_identifier, round_idx)
            catch e
                @warn "Could not read WT distance matrix for $(config.wt_identifier) round $round_idx: $e. Using identity matrix."
                dist_wt = Matrix{Float64}(I, size(PAE_wt))
            end
        end
        
        # Calculate terms for each index in this round
        for (term_idx, i) in enumerate(indices)
            if matrix_idx <= size(dist_mut, 1) && i <= size(dist_mut, 2) && 
               matrix_idx <= size(dist_wt, 1) && i <= size(dist_wt, 2) &&
               matrix_idx <= size(PAE_mut, 1) && i <= size(PAE_mut, 2) &&
               matrix_idx <= size(PAE_wt, 1) && i <= size(PAE_wt, 2)
                
                d_mut = dist_mut[matrix_idx, i]
                d_wt = dist_wt[matrix_idx, i]
                
                if d_mut > 0.0 && d_wt > 0.0
                    # Store individual terms for averaging
                    mut_terms[round_idx, term_idx] = PAE_mut[matrix_idx, i]^(2-rho) / (d_mut^α)
                    wt_terms[round_idx, term_idx] = PAE_wt[matrix_idx, i]^(2-rho) / (d_wt^α)
                end
            end
        end
    end
    
    # Average across rounds first, then subtract (as requested)
    avg_mut_terms = mean(mut_terms, dims=1)[1, :]  # Average over rounds
    avg_wt_terms = mean(wt_terms, dims=1)[1, :]    # Average over rounds
    
    # Calculate ΔΔS using the averaged terms
    ΔΔS = 0.0
    for (term_idx, i) in enumerate(indices)
        ΔΔS += abs(Γ[matrix_idx, i]) * (avg_mut_terms[term_idx] - avg_wt_terms[term_idx])
    end
    
    ΔΔS = ΔΔS / length(indices)
    return ΔΔS
end

function ΔΔG_prime(A::Float64, ΔΔS::Float64, ΔΔG::Float64)
    return ΔΔG + A * ΔΔS
end

"""
    calculate_ddgs(task_file_path::String, single_ddG::Dict{String, Float64}, pdb_path::String, WT_pae::Vector{Matrix{Float64}}, paes::Dict{String, Vector{Matrix{Float64}}}, ddG_exp::DataFrame, rho::Float64, A::Float64, config::EntropyConfig; verbose::Bool=false)

Calculate the predicted ΔΔG values for a set of mutations using multi-round data.

# Arguments
- `task_file_path::String`: Path to file containing list of mutations to analyze
- `single_ddG::Dict{String, Float64}`: Dictionary mapping mutation strings to Rosetta ddG values
- `pdb_path::String`: Path to the PDB file for computing the Gamma matrix
- `WT_pae::Vector{Matrix{Float64}}`: Vector of PAE matrices for wild-type protein (one per round)
- `paes::Dict{String, Vector{Matrix{Float64}}}`: Dictionary mapping mutations to their vectors of PAE matrices
- `ddG_exp::DataFrame`: DataFrame containing experimental ddG values
- `rho::Float64`: Parameter controlling contribution of PAE differences
- `A::Float64`: Scaling parameter for entropy contribution
- `α::Float64`: Exponent for distance scaling in entropy calculation
- `config::EntropyConfig`: Configuration for entropy calculations
- `verbose::Bool`: If true, prints details about skipped variants

# Returns
- `filtered_ddG_exp::Vector{Float64}`: Filtered experimental ddG values
- `ΔΔGs::Vector{Float64}`: Predicted total ddG values including entropy contribution
- `r_ddGs::Vector{Float64}`: Original Rosetta ddG predictions
"""
function calculate_ddgs(task_file_path::String, single_ddG::Dict{String, Float64}, pdb_path::String, WT_pae::Vector{Matrix{Float64}}, paes::Dict{String, Vector{Matrix{Float64}}}, ddG_exp::DataFrame, rho::Float64, A::Float64, α::Float64, config::EntropyConfig; verbose::Bool=false)
    # Compute Gamma matrix from PDB file
    coordinates = read_coordinates(pdb_path)
    Γ = compute_Γ(coordinates)
    
    # Process mutations
    mutations = read_mutations_from_file(task_file_path)
    results = Vector{Tuple{Float64, Float64, Float64}}()
    
    for m in mutations
        position = parse_mutation_position(m)
        result = process_single_mutation(m, position, single_ddG, paes, Γ, WT_pae, ddG_exp, rho, A, α, config, verbose)
        if result !== nothing && result[2] !== NaN
            push!(results, result)
        else
            if verbose == true
                println("Skipping mutation: $m (result is NaN or not found)")
            end
        end
    end
    
    # Extract results
    isempty(results) && return Float64[], Float64[], Float64[]
    return [r[1] for r in results], [r[2] for r in results], [r[3] for r in results]
end

"""Extract the position number from a mutation string."""
parse_mutation_position(mutation::AbstractString) = parse(Int, match(r"\d+", mutation).match)

"""Get the experimental ΔΔG value for a specific position and mutation."""
function get_experimental_ddg(ddG_exp::DataFrame, position::Int, mutation_residue::String)
    matching_rows = (ddG_exp.position .== position) .& (ddG_exp.mutation .== mutation_residue)
    return mean(ddG_exp[matching_rows, :ddG])
end

"""Read mutation list from file and return as vector."""
function read_mutations_from_file(task_file_path::String)
    return [strip(line) for line in eachline(task_file_path)]
end

"""Process a single mutation and return the calculated values."""
function process_single_mutation(mutation::AbstractString, position::Int, single_ddG::Dict{String, Float64}, 
                                paes::Dict{String, Vector{Matrix{Float64}}}, Γ::Matrix{Float64}, 
                                WT_pae::Vector{Matrix{Float64}}, ddG_exp::DataFrame, rho::Float64, A::Float64, α::Float64, 
                                config::EntropyConfig, verbose::Bool=false)
    
    mutation_upper = uppercase(mutation)
    if !haskey(single_ddG, mutation_upper)
        if verbose
            println("Skipping variant: $mutation (not found in single_ddG)")
        end
        return nothing
    end
    
    if !haskey(paes, mutation)
        if verbose
            println("Skipping variant: $mutation (no PAE data found)")
        end
        return nothing
    end
    
    ddG = single_ddG[mutation_upper]
    pae_rounds = paes[mutation]
    ΔΔS_val = ΔΔS(position, rho, α, Γ, pae_rounds, WT_pae, mutation, config)
    predicted_ddG = ΔΔG_prime(A, ΔΔS_val, ddG)
    experimental_ddG = get_experimental_ddg(ddG_exp, position, string(mutation_upper[end]))
    
    return (experimental_ddG, predicted_ddG, ddG)
end