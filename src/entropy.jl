"""
    find_stable_neighbors(matrix_idx::Int, distance_matrices::Vector{Matrix{Float64}}, mutation::String)

Find residues that are within 13Å in all rounds using matrix indices.

Args:
- matrix_idx: Position index in the truncated matrix (1-based)
- distance_matrices: Vector of truncated distance matrices  
- mutation: Mutation identifier for logging

Returns:
- Vector of matrix indices for stable neighbor residues
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
                    # Convert PDB residue numbers to matrix indices
                    # PDB residue number -> matrix index using the offset from params
                    # matrix_index = pdb_resnum - offset
                    first_pdb_in_matrix = params.offset + 1  # First PDB residue in truncated matrix
                    low_plddt_matrix_indices = [r - params.offset for r in low_plddt_residues if r >= first_pdb_in_matrix]
                    union!(all_low_plddt, low_plddt_matrix_indices)
                end
            catch e
                # Skip rounds that don't exist or have errors - this is expected behavior
                continue
            end
        end
    end

    filtered_indices = filter(idx -> !(idx in all_low_plddt), indices)
    
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

The function converts biological residue numbers to matrix indices and performs entropy calculations
on truncated matrices where the first matrix position corresponds to the first biological residue
in the analysis region.
"""
function ΔΔS(params::EntropyParams, data::MutationData, self_normalize::Bool)::Float64
    # Convert biological residue number to matrix index for truncated matrices
    # For Thermonuclease: position=88, offset=87 → matrix_idx=1 (first position in truncated matrix)
    matrix_idx = params.position - params.offset
    
    # First layer filtering: check if the current position itself has low pLDDT
    if params.filter_low_plddt
        if !isempty(params.data_dir)
            # Use the biological residue number for pLDDT lookup
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
    
    # Find stable neighbors for both mutant and WT using matrix indices
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
        @warn "All residues within 13Å were filtered out due to low pLDDT for position $(params.position) in $(data.mutation)."
        return NaN
    end
    
    # Calculate entropy terms
    mut_terms, wt_terms = calculate_entropy_terms(matrix_idx, indices, round_data, params)
    
    # Calculate final result
    avg_mut_terms = mean(mut_terms, dims = 1)[1, :]
    avg_wt_terms = mean(wt_terms, dims = 1)[1, :]
    
    ΔΔS_val = sum(avg_mut_terms - avg_wt_terms)
    if self_normalize
        ΔΔS_val /= length(indices)
    end
    return ΔΔS_val
end

"""
    ΔΔG_prime(A, ΔΔS, ΔΔG)

Calculates the corrected ΔΔG value.
"""
ΔΔG_prime(A::Float64, ΔΔS::Float64, ΔΔG::Float64)::Float64 = ΔΔG + A * ΔΔS


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
    self_normalize::Bool = false,
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
    
    ΔΔS_val = ΔΔS(params, data, self_normalize)
    predicted_ddG = ΔΔG_prime(A, ΔΔS_val, ddG)
    
    return (experimental_ddG, predicted_ddG, ddG)
end

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
    self_normalize::Bool = false,
)::Tuple{Vector{Float64},Vector{Float64},Vector{Float64}}
    mutations = read_mutations_from_file(task_file_path)
    results = Vector{Tuple{Float64,Float64,Float64}}()

    for m in mutations
        position = parse_mutation_position(m)
        result = process_single_mutation(m, position, single_ddG, wt_pae, wt_dist, paes, dist_matrices, ddG_exp, rho, A, α, offset, verbose, data_dir, filter_low_plddt, plddt_threshold, self_normalize)
        
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

"""
    process_entropy_data(datadir::String, param_subdir::String, nearby_normalize::Bool, verbose::Bool=false)::DataFrame

Process entropy data from saved JLD2 files.

Reads JLD2 files with the following structure:
- metadata["mutation_position"]: biological residue number of mutation
- metadata["indexing_info"]["matrix_offset"]: offset for converting to matrix index
- results["all_residues"]["ddS_filtered"]: ddS values array
- results["nearby_residues"]["ddS_filtered"]: ddS values for nearby residues

Matrix index is calculated as: mutation_position - matrix_offset

# Arguments
- `datadir`: Base directory containing mutant subdirectories
- `param_subdir`: Parameter subdirectory name (e.g., "alpha_0.5_rho_1.0")
- `nearby_normalize`: If true, divide nearby_ddS by count of nearby residues
- `verbose`: If true, print progress and warning messages

# Returns
- DataFrame with columns: mutant, mutant_ddS, nearby_ddS
"""
function process_entropy_data(datadir::String, param_subdir::String, nearby_normalize::Bool, verbose::Bool=false)::DataFrame
    all_mutant_data = Dict{String, Dict{String, Float64}}()

    # Iterate through each mutant directory in the base data directory
    for mutant_dir_name in readdir(datadir)
        mutant_path = joinpath(datadir, mutant_dir_name)
        if !isdir(mutant_path)
            continue
        end

        # Construct the path to the specific parameter subdirectory
        param_path = joinpath(mutant_path, param_subdir)
        
        if !isdir(param_path)
            if verbose
                println("Parameter directory not found for mutant $mutant_dir_name: $param_path. Skipping...")
            end
            continue
        end

        # Find the JLD2 file in the parameter directory
        jld_files = filter(f -> endswith(f, ".jld2"), readdir(param_path))
        if isempty(jld_files)
            if verbose
                println("No JLD2 file found in $param_path. Skipping...")
            end
            continue
        end
        filepath = joinpath(param_path, first(jld_files))

        # Load data from the JLD2 file
        local temp
        try
            temp = jldopen(filepath, "r") do file
                if !haskey(file, "data")
                    println("JLD2 file $filepath is missing the 'data' key. Skipping...")
                    return nothing
                end
                file["data"]
            end
            if temp === nothing
                continue
            end
        catch e
            println("Failed to open or read JLD2 file $filepath. Error: $e. Skipping...")
            continue
        end

        # Extract metadata
        meta = temp["metadata"]
        mutant_name = uppercase(meta["mutant"])  # Use uppercase for consistency with experimental data
        
        # Validate data format - must have indexing_info
        if !haskey(meta, "indexing_info")
            println("Invalid data format for mutant $mutant_name: missing 'indexing_info'. Skipping...")
            continue
        end
        
        # Get mutation position and convert to matrix index
        mutation_position = meta["mutation_position"]
        matrix_offset = meta["indexing_info"]["matrix_offset"]
        mutation_matrix_idx = mutation_position - matrix_offset
        
        if verbose
            println("Processing $mutant_name: position=$mutation_position, matrix_offset=$matrix_offset, matrix_idx=$mutation_matrix_idx")
        end

        # Extract results data
        results = temp["results"]

        # Calculate the ddS at mutation site
        local mutant_ddS = 0.0
        ddS_array = results["all_residues"]["ddS_filtered"]
        if mutation_matrix_idx >= 1 && mutation_matrix_idx <= length(ddS_array)
            val = ddS_array[mutation_matrix_idx]
            if !isnan(val)
                mutant_ddS = val
            elseif verbose
                println("Warning: NaN ddS value at mutation site for $mutant_name")
            end
        else
            if verbose
                println("Warning: Matrix index $mutation_matrix_idx out of bounds for $mutant_name (array length: $(length(ddS_array)))")
            end
        end

        # Calculate the ddS sum for nearby residues
        local nearby_ddS = 0.0
        local nearby_count = 0
        
        nearby_residues_data = results["nearby_residues"]
        if haskey(nearby_residues_data, "ddS_filtered")
            nearby_ddS_values = filter(!isnan, nearby_residues_data["ddS_filtered"])
            nearby_ddS = sum(nearby_ddS_values)
            nearby_count = length(nearby_ddS_values)
        end
        
        if verbose && nearby_count == 0
            println("Warning: No valid nearby residues found for mutant $mutant_name")
        end
        
        # Normalize by count if flag is set
        if nearby_normalize && nearby_count > 0
            nearby_ddS = nearby_ddS / nearby_count
        end

        # Store the results in the dictionary
        all_mutant_data[mutant_name] = Dict(
            "mutant_ddS" => mutant_ddS,
            "nearby_ddS" => nearby_ddS
        )
    end

    # Convert the calculated data dictionary to a DataFrame
    calculated_df = DataFrame(
        mutant = collect(keys(all_mutant_data)),
        mutant_ddS = [d["mutant_ddS"] for d in values(all_mutant_data)],
        nearby_ddS = [d["nearby_ddS"] for d in values(all_mutant_data)]
    )

    return calculated_df
end