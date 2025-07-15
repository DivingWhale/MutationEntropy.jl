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