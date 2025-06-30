"""
    Φ(r, η, κ)

Computes the pairwise interaction strength based on a generalized exponential function.
"""
Φ(r::Real, η::Real, κ::Real)::Float64 = exp(-(r / η)^κ)

"""
    Γ(coordinates, η, κ)

Constructs the Γ matrix, a key component in modeling protein dynamics, from atomic coordinates.

# Arguments
- `coordinates::AbstractVector{<:AbstractVector{<:Real}}`: A vector of 3D coordinates for each atom.
- `η::Integer`: Scaling parameter for the distance.
- `κ::Integer`: Exponent parameter determining the shape of the interaction falloff.

# Returns
- `Matrix{Float64}`: The computed Γ matrix.
"""
function Γ(coordinates::AbstractVector{<:AbstractVector{<:Real}}, η::Integer, κ::Integer)::Matrix{Float64}
    natoms = length(coordinates)
    Γ_matrix = zeros(Float64, natoms, natoms)

    for i in 1:natoms, j in (i+1):natoms
        dist = norm(coordinates[i] - coordinates[j])
        interaction = -Φ(dist, η, κ)
        Γ_matrix[i, j] = interaction
        Γ_matrix[j, i] = interaction  # Exploiting symmetry
    end

    # Diagonal elements are the negative sum of the off-diagonal elements in each column
    for i in 1:natoms
        Γ_matrix[i, i] = -sum(view(Γ_matrix, :, i))
    end

    return Γ_matrix
end

"""
    compute_Γ(coordinates)

Computes the final Γ matrix by combining results from two different parameter sets.
"""
function compute_Γ(coordinates::AbstractVector{<:AbstractVector{<:Real}})::Matrix{Float64}
    Γ1 = Γ(coordinates, 20, 7)
    Γ2 = Γ(coordinates, 13, 10)
    return Γ1 + Γ2
end

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
    msf(Γ)

Calculates the mean square fluctuation from the Γ matrix.
"""
msf(Γ::AbstractMatrix{<:Real})::Vector{Float64} = diag(pinv(Γ))

"""
    read_single_pae(data_path, mutation, round_val)

Reads a single Predicted Aligned Error (PAE) matrix from a specified directory.
"""
function read_single_pae(data_path::String, mutation::AbstractString, round_val::Int)::Union{Matrix{Float64}, Nothing}
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

Reads all PAE matrices for a list of mutations, utilizing a cache to avoid redundant reads.
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
                round_pae = read_single_pae(source_data_path, mutation, round_idx)
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
    EntropyConfig

Configuration for entropy calculations.
"""
struct EntropyConfig
    datadir::String
    wt_identifier::String
    offset::Int

    function EntropyConfig(
        datadir::String;
        wt_identifier::String = "WT",
        offset::Int = 0,
    )
        new(datadir, wt_identifier, offset)
    end
end

"""
    ΔΔS(position, rho, α, Γ, PAE_mut_rounds, PAE_wt_rounds, mutation, config; ...)

Calculates the change in entropy (ΔΔS) for a given mutation.
"""
function ΔΔS(
    position::Int,
    rho::Float64,
    Γ::Matrix{Float64},
    PAE_mut_rounds::Vector{Matrix{Float64}},
    PAE_wt_rounds::Vector{Matrix{Float64}},
    mutation::AbstractString,
    config::EntropyConfig;
    cache_manager::Union{CacheManager,Nothing} = nothing,
)::Float64
    matrix_idx = position - config.offset
    indices = findall(x -> abs(x) > 1e-5, Γ[matrix_idx, :])
    num_rounds = min(length(PAE_mut_rounds), length(PAE_wt_rounds))

    if num_rounds == 0
        @warn "No matching rounds for mutation $mutation."
        return 0.0
    end

    mut_terms = zeros(Float64, num_rounds, length(indices))
    wt_terms = zeros(Float64, num_rounds, length(indices))

    for round_idx in 1:num_rounds
        PAE_mut = PAE_mut_rounds[round_idx]
        PAE_wt = PAE_wt_rounds[round_idx]

        for (term_idx, i) in enumerate(indices)
            if checkbounds(Bool, PAE_mut, matrix_idx, i) &&
               checkbounds(Bool, PAE_wt, matrix_idx, i)
                
                mut_terms[round_idx, term_idx] = PAE_mut[matrix_idx, i]^(2 - rho)
                wt_terms[round_idx, term_idx] = PAE_wt[matrix_idx, i]^(2 - rho)
            end
        end
    end

    avg_mut_terms = mean(mut_terms, dims = 1)[1, :]
    avg_wt_terms = mean(wt_terms, dims = 1)[1, :]

    ΔΔS_val = 0.0
    for (term_idx, i) in enumerate(indices)
        ΔΔS_val += abs(Γ[matrix_idx, i]) * (avg_mut_terms[term_idx] - avg_wt_terms[term_idx])
    end

    return ΔΔS_val / length(indices)
end

"""
    ΔΔG_prime(A, ΔΔS, ΔΔG)

Calculates the corrected ΔΔG value.
"""
ΔΔG_prime(A::Float64, ΔΔS::Float64, ΔΔG::Float64)::Float64 = ΔΔG + A * ΔΔS

"""
    calculate_ddgs(...)

Calculates predicted ΔΔG values for a set of mutations.
"""
function calculate_ddgs(
    task_file_path::String,
    single_ddG::Dict{String,Float64},
    pdb_path::String,
    WT_pae::Vector{Matrix{Float64}},
    paes::Dict{String,Vector{Matrix{Float64}}},
    ddG_exp::DataFrame,
    rho::Float64,
    A::Float64,
    config::EntropyConfig;
    verbose::Bool = false,
    cache_manager::Union{CacheManager,Nothing} = nothing,
)::Tuple{Vector{Float64},Vector{Float64},Vector{Float64}}
    coordinates = read_coordinates(pdb_path)
    Γ = compute_Γ(coordinates)
    mutations = read_mutations_from_file(task_file_path)
    results = Vector{Tuple{Float64,Float64,Float64}}()

    for m in mutations
        position = parse_mutation_position(m)
        result = process_single_mutation(m, position, single_ddG, paes, Γ, WT_pae, ddG_exp, rho, A, config, verbose; cache_manager)
        
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

"""Processes a single mutation to calculate its ΔΔG."""
function process_single_mutation(
    mutation::AbstractString,
    position::Int,
    single_ddG::Dict{String,Float64},
    paes::Dict{String,Vector{Matrix{Float64}}},
    Γ::Matrix{Float64},
    WT_pae::Vector{Matrix{Float64}},
    ddG_exp::DataFrame,
    rho::Float64,
    A::Float64,
    config::EntropyConfig,
    verbose::Bool;
    cache_manager::Union{CacheManager,Nothing},
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

    experimental_ddG = get_experimental_ddg(ddG_exp, position, string(mutation_upper[end]))
    if experimental_ddG === nothing
        verbose && println("Skipping: $mutation (no experimental ddG).")
        return nothing
    end

    ddG = single_ddG[mutation_upper]
    pae_rounds = paes[mutation]
    
    ΔΔS_val = ΔΔS(position, rho, Γ, pae_rounds, WT_pae, mutation, config; cache_manager)
    predicted_ddG = ΔΔG_prime(A, ΔΔS_val, ddG)
    
    return (experimental_ddG, predicted_ddG, ddG)
end
