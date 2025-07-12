"""
    read_coordinates(filename)

Reads CA atom coordinates from a PDB file.
"""
function read_coordinates(filename::String)::Vector{Vector{Float64}}
    structure = read(filename, PDBFormat)
    ca_atoms = collectatoms(structure, sel"name CA")
    return getfield.(ca_atoms, :coords)
end

"""
    read_xvg(filename)

Parses a generic .xvg file, returning a vector of (Int, Float64) tuples.
"""
function read_xvg(filename::String)::Vector{Tuple{Int, Float64}}
    results = Tuple{Int, Float64}[]
    for line in eachline(filename)
        if startswith(line, ('#', '@')) || isempty(line)
            continue
        end
        parts = split(line)
        if length(parts) >= 2
            push!(results, (parse(Int, parts[1]), parse(Float64, parts[2])))
        end
    end
    return results
end



"""
    read_round_pae(data_path, mutation, round_val)

Reads a single Predicted Aligned Error (PAE) matrix for a specific round from a directory.
"""
function read_round_pae(data_path::String, mutation::AbstractString, round_val::Int)::Union{Matrix{Float64}, Nothing}
    round_dir = joinpath(data_path, "$(mutation)_$(round_val)")
    !isdir(round_dir) && return nothing

    seed_dirs = sort(glob("seed-*_sample-0", round_dir))
    isempty(seed_dirs) && return nothing

    confidence_file = joinpath(first(seed_dirs), "confidences.json")
    !isfile(confidence_file) && return nothing

    try
        raw_data = JSON.parsefile(confidence_file)
        pae_data = raw_data["pae"]
        # Convert vector of vectors to a matrix and ensure Float64 type
        return reduce(hcat, [map(Float64, row) for row in pae_data])'
    catch e
        @warn "Failed to read or parse PAE data from $confidence_file: $e"
        return nothing
    end
end

"""
    read_paes(task_file_path, source_data_path; ...)

Reads all PAE matrices for a list of mutations from a task file, utilizing a cache.
"""
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

    # Check cache for existing data
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

    # Process uncached mutations in parallel
    if !isempty(mutations_to_process)
        results_channel = Channel{Tuple{String,Vector{Matrix{Float64}}}}(length(mutations_to_process))

        Threads.@threads for mutation in mutations_to_process
            pae_matrices = Vector{Matrix{Float64}}()
            for round_idx in 1:num_rounds
                round_pae = read_round_pae(source_data_path, mutation, round_idx)
                round_pae !== nothing && push!(pae_matrices, round_pae)
            end

            if !isempty(pae_matrices)
                put!(results_channel, (mutation, pae_matrices))
            elseif verbose
                @warn "No PAE data found for $mutation."
            end
        end
        close(results_channel)

        # Collect results and update cache
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

"""
    get_variant_paes(mutation, source_data_path; ...)

Gets all PAE matrices for a single variant, utilizing a cache.
"""
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

    # Check cache for existing data
    cached_pae = load_cached_data(cache_manager, mutation, :pae)
    if cached_pae !== nothing
        verbose && println("Loaded $mutation from cache.")
        return cached_pae, cache_manager
    end

    # If not in cache, read from source
    verbose && println("Cache miss for $mutation. Reading from source.")
    pae_matrices = Vector{Matrix{Float64}}()
    for round_idx in 1:num_rounds
        round_pae = read_round_pae(source_data_path, mutation, round_idx)
        round_pae !== nothing && push!(pae_matrices, round_pae)
    end

    if isempty(pae_matrices)
        verbose && @warn "No PAE data found for $mutation."
        return nothing
    end

    # Save to cache and return
    verbose && println("Saving $mutation to cache.")
    save_cached_data(cache_manager, mutation, :pae, pae_matrices)
    
    return pae_matrices, cache_manager
end

"""
    MutationData

Contains all matrices and identifiers for a single mutation analysis.
"""
struct MutationData
    wt_pae::Vector{Matrix{Float64}}
    mutant_pae::Vector{Matrix{Float64}}
    wt_dist::Vector{Matrix{Float64}}
    mutant_dist::Vector{Matrix{Float64}}
    mutation::String
end

"""
    EntropyParams

Parameters for entropy calculation.
"""
struct EntropyParams
    position::Int
    rho::Float64
    α::Float64
    offset::Int
    filter_low_plddt::Bool
    plddt_threshold::Float64
    data_dir::String
end



"""
    find_stable_neighbors(matrix_idx::Int, distance_matrices::Vector{Matrix{Float64}}, mutation::String)

Find residues that are within 13Å in all rounds.
"""
function find_stable_neighbors(matrix_idx::Int, distance_matrices::Vector{Matrix{Float64}}, mutation::String)::Vector{Int}
    all_round_indices = Vector{Vector{Int}}()
    
    for dist_matrix in distance_matrices
        round_indices = find_residues_within_distance(matrix_idx, dist_matrix; distance=13.0)
        push!(all_round_indices, round_indices)
    end
    
    if isempty(all_round_indices)
        @warn "No valid distance matrices found for $mutation."
        return Int[]
    end
    
    # Get intersection of all rounds
    indices = all_round_indices[1]
    for round_indices in all_round_indices[2:end]
        indices = intersect(indices, round_indices)
    end
    
    return indices
end

"""
    filter_low_plddt_residues_per_round(indices::Vector{Int}, mutation::String, params::EntropyParams, num_rounds::Int=20)

Filter out low pLDDT residues for each round and return the union of all filtered residues.
Uses data_dir from params struct.
"""
function filter_low_plddt_residues_per_round(indices::Vector{Int}, mutation::String, params::EntropyParams, num_rounds::Int=20)::Vector{Int}
    if !params.filter_low_plddt
        return indices
    end
    
    # Ensure data_dir is provided
    if isempty(params.data_dir)
        error("No data directory provided: params.data_dir is empty. Please provide data_dir in EntropyParams.")
    end
    
    all_low_plddt = Set{Int}()
    datadir = params.data_dir
    
    # Get low pLDDT residues for both mutant and WT for each available round
    for mut_id in [mutation, "WT"]
        for round_idx in 1:num_rounds
            try
                low_plddt_residues = get_low_plddt_residues(mut_id, round_idx, datadir, threshold=params.plddt_threshold)
                if !isempty(low_plddt_residues)
                    # Convert full sequence residue numbers to matrix indices
                    # Full sequence residue 88 -> matrix index 1, 89 -> 2, etc.
                    # So matrix_index = full_seq_residue - 87
                    low_plddt_matrix_indices = [r - 87 for r in low_plddt_residues if r >= 88]  # Only include residues in PDB region
                    union!(all_low_plddt, low_plddt_matrix_indices)
                end
            catch e
                # Skip rounds that don't exist or have errors - this is expected behavior
                continue
            end
        end
    end
    
    # Filter out low pLDDT residues from indices
    # original_count = length(indices)
    filtered_indices = filter(idx -> !(idx in all_low_plddt), indices)
    # filtered_count = original_count - length(filtered_indices)
    
    # if filtered_count > 0
    #     println("ΔΔS calculation: Filtered out $filtered_count low pLDDT residues (pLDDT < $(params.plddt_threshold)) across all rounds for position $(params.position)")
    # end
    
    return filtered_indices
end

"""
    calculate_entropy_terms(matrix_idx::Int, indices::Vector{Int}, round_data::Vector{Tuple{Matrix{Float64}, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}}}, params::EntropyParams)

Calculate entropy terms for mutant and wild-type across all rounds.
"""
function calculate_entropy_terms(matrix_idx::Int, indices::Vector{Int}, round_data::Vector{Tuple{Matrix{Float64}, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}}}, params::EntropyParams)
    mut_terms = zeros(length(round_data), length(indices))
    wt_terms = zeros(length(round_data), length(indices))
    
    for (round_idx, (PAE_mut, PAE_wt, dist_mut, dist_wt)) in enumerate(round_data)
        for (term_idx, i) in enumerate(indices)
            if checkbounds(Bool, dist_mut, matrix_idx, i) &&
               checkbounds(Bool, dist_wt, matrix_idx, i) &&
               checkbounds(Bool, PAE_mut, matrix_idx, i) &&
               checkbounds(Bool, PAE_wt, matrix_idx, i)
                
                d_mut = dist_mut[matrix_idx, i]
                d_wt = dist_wt[matrix_idx, i]

                if d_mut > 0.0 && d_wt > 0.0
                    mut_terms[round_idx, term_idx] = PAE_mut[matrix_idx, i]^(2 - params.rho) / (d_mut^params.α)
                    wt_terms[round_idx, term_idx] = PAE_wt[matrix_idx, i]^(2 - params.rho) / (d_wt^params.α)
                end
            end
        end
    end
    
    return mut_terms, wt_terms
end

"""
    ΔΔS(params::EntropyParams, data::MutationData)

Calculates the change in entropy (ΔΔS) for a given mutation using structured parameters.
This is the main function that orchestrates the entropy calculation process.
"""
function ΔΔS(params::EntropyParams, data::MutationData)::Float64
    matrix_idx = params.position - params.offset
    
    # First layer filtering: check if the current position itself is low pLDDT in either WT or mutant
    if params.filter_low_plddt
        if !isempty(params.data_dir)
            full_seq_residue_number = params.position
            is_position_low_plddt = false
            
            # Check both mutant and WT pLDDTs
            for mut_id in [data.mutation, "WT"]
                num_rounds_to_check = (mut_id == "WT") ? length(data.wt_pae) : length(data.mutant_pae)
                
                for round_idx in 1:num_rounds_to_check
                    try
                        low_plddt_residues = get_low_plddt_residues(mut_id, round_idx, params.data_dir, threshold=params.plddt_threshold)
                        if full_seq_residue_number in low_plddt_residues
                            is_position_low_plddt = true
                            break
                        end
                    catch e
                        # Skip rounds that don't exist - this is expected
                        continue
                    end
                end
                if is_position_low_plddt
                    break
                end
            end
            
            if is_position_low_plddt
                # println("ΔΔS calculation: Position $(params.position) itself is low pLDDT (< $(params.plddt_threshold)) in $(data.mutation) or WT")
                return NaN
            end
        end
    end
    
    # Collect round data
    round_data = Vector{Tuple{Matrix{Float64}, Matrix{Float64}, Matrix{Float64}, Matrix{Float64}}}()
    for round_idx in eachindex(data.mutant_pae, data.wt_pae)
        if round_idx <= length(data.mutant_dist) && round_idx <= length(data.wt_dist)
            push!(round_data, (data.mutant_pae[round_idx], data.wt_pae[round_idx], 
                              data.mutant_dist[round_idx], data.wt_dist[round_idx]))
        end
    end
    
    if isempty(round_data)
        @warn "No valid round data found for position $(params.position) in $(data.mutation)."
        return NaN
    end
    
    # Find stable neighbors for both mutant and WT
    indices_mut = find_stable_neighbors(matrix_idx, data.mutant_dist, "$(data.mutation)_mutant")
    indices_wt = find_stable_neighbors(matrix_idx, data.wt_dist, "$(data.mutation)_wt")
    
    # Use intersection of both mutant and WT stable neighbors
    indices = intersect(indices_mut, indices_wt)
    
    if isempty(indices)
        @warn "No residues found within 13Å in all rounds for position $(params.position) in $(data.mutation)."
        return NaN
    end
    
    # Filter out low pLDDT residues if requested (applies to all rounds)
    indices = filter_low_plddt_residues_per_round(indices, data.mutation, params, length(round_data))
    
    if isempty(indices)
        # @warn "All residues within 13Å were filtered out due to low pLDDT for position $(params.position) in $(data.mutation)."
        return NaN
    end
    
    # Calculate entropy terms
    mut_terms, wt_terms = calculate_entropy_terms(matrix_idx, indices, round_data, params)
    
    # Calculate final result
    avg_mut_terms = mean(mut_terms, dims = 1)[1, :]
    avg_wt_terms = mean(wt_terms, dims = 1)[1, :]
    
    ΔΔS_val = sum(avg_mut_terms - avg_wt_terms)
    return ΔΔS_val
end

"""
    ΔΔG_prime(A, ΔΔS, ΔΔG)

Calculates the corrected ΔΔG value.
"""
ΔΔG_prime(A::Float64, ΔΔS::Float64, ΔΔG::Float64)::Float64 = ΔΔG + A * ΔΔS

"""
    calculate_ddgs(task_file_path, single_ddG, pdb_path, wt_pae, wt_dist, paes, dist_matrices, ddG_exp, rho, A, α, offset; ...)

Calculates predicted ΔΔG values for a set of mutations using the new structured API.
"""
function calculate_ddgs(
    task_file_path::String,
    single_ddG::Dict{String,Float64},
    pdb_path::String,
    wt_pae::Vector{Matrix{Float64}},
    wt_dist::Vector{Matrix{Float64}},
    paes::Dict{String,Vector{Matrix{Float64}}},
    dist_matrices::Dict{String,Vector{Matrix{Float64}}},
    ddG_exp::DataFrame,
    rho::Float64,
    A::Float64,
    α::Float64,
    offset::Int;
    verbose::Bool = false,
    data_dir::String = "",
    filter_low_plddt::Bool = false,
    plddt_threshold::Float64 = 90.0,
)::Tuple{Vector{Float64},Vector{Float64},Vector{Float64}}
    mutations = read_mutations_from_file(task_file_path)
    results = Vector{Tuple{Float64,Float64,Float64}}()

    for m in mutations
        position = parse_mutation_position(m)
        result = process_single_mutation(m, position, single_ddG, wt_pae, wt_dist, paes, dist_matrices, ddG_exp, rho, A, α, offset, verbose, data_dir, filter_low_plddt, plddt_threshold)
        
        if result !== nothing && !isnan(last(result))
            push!(results, result)
        elseif verbose
            println("Skipping mutation: $m (result is NaN or not found).")
        end
    end

    isempty(results) && return (Float64[], Float64[], Float64[])
    
    # Unzip results into separate vectors
    exp_ddG = [r[1] for r in results]
    pred_ddG = [r[2] for r in results]
    rosetta_ddG = [r[3] for r in results]
    
    return exp_ddG, pred_ddG, rosetta_ddG
end

"""Extracts the position number from a mutation string (e.g., "A123G" -> 123)."""
parse_mutation_position(mutation::AbstractString)::Int = parse(Int, match(r"\d+", mutation).match)

"""Retrieves the experimental ΔΔG value for a mutation."""
function get_experimental_ddg(ddG_exp::DataFrame, position::Int, mutation_residue::String)::Union{Float64, Nothing}
    matching_rows = (ddG_exp.position .== position) .& (ddG_exp.mutation .== mutation_residue)
    
    if !any(matching_rows)
        return nothing
    end
    
    return mean(ddG_exp[matching_rows, :ddG])
end


"""Reads a list of mutations from a file."""
function read_mutations_from_file(task_file_path::String)::Vector{String}
    return [strip(line) for line in eachline(task_file_path) if !isempty(strip(line))]
end

"""Processes a single mutation to calculate its ΔΔG using structured data."""
function process_single_mutation(
    mutation::AbstractString,
    position::Int,
    single_ddG::Dict{String,Float64},
    wt_pae::Vector{Matrix{Float64}},
    wt_dist::Vector{Matrix{Float64}},
    paes::Dict{String,Vector{Matrix{Float64}}},
    dist_matrices::Dict{String,Vector{Matrix{Float64}}},
    ddG_exp::DataFrame,
    rho::Float64,
    A::Float64,
    α::Float64,
    offset::Int,
    verbose::Bool,
    data_dir::String = "",
    filter_low_plddt::Bool = false,
    plddt_threshold::Float64 = 90.0,
)::Union{Tuple{Float64,Float64,Float64},Nothing}
    mutation_upper = uppercase(mutation)
    
    if !haskey(single_ddG, mutation_upper)
        verbose && println("Skipping: $mutation (not in Rosetta ddG data).")
        return nothing
    end
    
    if !haskey(paes, mutation)
        verbose && println("Skipping: $mutation (no PAE data).")
        return nothing
    end
    
    if !haskey(dist_matrices, mutation)
        verbose && println("Skipping: $mutation (no distance matrix data).")
        return nothing
    end

    experimental_ddG = get_experimental_ddg(ddG_exp, position, string(mutation_upper[end]))
    if experimental_ddG === nothing
        verbose && println("Skipping: $mutation (no experimental ddG).")
        return nothing
    end

    ddG = single_ddG[mutation_upper]
    
    # Create structured data
    data = MutationData(wt_pae, paes[mutation], wt_dist, dist_matrices[mutation], mutation)
    params = EntropyParams(position, rho, α, offset, filter_low_plddt, plddt_threshold, data_dir)
    
    ΔΔS_val = ΔΔS(params, data)
    predicted_ddG = ΔΔG_prime(A, ΔΔS_val, ddG)
    
    return (experimental_ddG, predicted_ddG, ddG)
end

"""
    get_residue_calpha_b_factor(pdb_file_path)

Reads C-alpha B-factors from a PDB file and returns a dictionary mapping residue identifiers to B-factor values.
"""
function get_residue_calpha_b_factor(pdb_file_path::String)::Dict{String,Float64}
    residue_calpha_b_factors = Dict{String,Float64}()
    try
        struc = read(pdb_file_path, PDBFormat)
        for residue in collectresidues(struc, standardselector)
            calpha_atom = nothing
            for atom_pair in atoms(residue)
                atom = atom_pair[2]  # Extract the atom from the pair
                if strip(String(atom_pair[1])) == "CA"  # atom_pair[1] is the atom name
                    calpha_atom = atom
                    break
                end
            end
            if calpha_atom !== nothing
                residue_id = "$(resname(residue)) $(resnumber(residue)) $(chainid(residue))"
                residue_calpha_b_factors[residue_id] = tempfactor(calpha_atom)
            end
        end
    catch e
        @warn "Error reading PDB file: $e"
    end
    return residue_calpha_b_factors
end