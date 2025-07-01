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
    compute_Γ(coordinates, config)

Computes the final Γ matrix by combining results from two different parameter sets using config parameters.
"""
function compute_Γ(coordinates::AbstractVector{<:AbstractVector{<:Real}}, config::EntropyConfig)::Matrix{Float64}
    Γ1 = Γ(coordinates, config.eta1, config.kappa1)
    Γ2 = Γ(coordinates, config.eta2, config.kappa2)
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
    wt_dist_matrices::Union{Vector{Matrix{Float64}},Nothing}
    mut_dist_matrices::Union{Dict{String,Vector{Matrix{Float64}}},Nothing}
    offset::Int
    eta1::Int
    kappa1::Int
    eta2::Int
    kappa2::Int

    function EntropyConfig(
        datadir::String,
        eta1::Int,
        kappa1::Int,
        eta2::Int,
        kappa2::Int;
        wt_identifier::String = "WT",
        wt_dist_matrices::Union{Vector{Matrix{Float64}},Nothing} = nothing,
        mut_dist_matrices::Union{Dict{String,Vector{Matrix{Float64}}},Nothing} = nothing,
        offset::Int = 0,
    )
        new(datadir, wt_identifier, wt_dist_matrices, mut_dist_matrices, offset, eta1, kappa1, eta2, kappa2)
    end
end

"""
    ΔΔS(position, rho, α, Γ, PAE_mut_rounds, PAE_wt_rounds, mutation, config; ...)

Calculates the change in entropy (ΔΔS) for a given mutation.
"""
function ΔΔS(
    position::Int,
    rho::Float64,
    α::Float64,
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

        dist_mut = if config.mut_dist_matrices !== nothing && 
                      haskey(config.mut_dist_matrices, String(mutation)) && 
                      round_idx <= length(config.mut_dist_matrices[String(mutation)])
            config.mut_dist_matrices[String(mutation)][round_idx]
        else
            try
                get_dist_map(config.datadir, String(mutation), round_idx; cache_manager)
            catch e
                @warn "Could not read mutant distance matrix for $mutation round $round_idx: $e. Using identity matrix."
                Matrix{Float64}(I, size(PAE_mut))
            end
        end

        dist_wt = if config.wt_dist_matrices !== nothing && round_idx <= length(config.wt_dist_matrices)
            config.wt_dist_matrices[round_idx]
        else
            try
                get_dist_map(config.datadir, config.wt_identifier, round_idx; cache_manager)
            catch e
                @warn "Could not read WT distance matrix for $(config.wt_identifier) round $round_idx: $e. Using identity matrix."
                Matrix{Float64}(I, size(PAE_wt))
            end
        end

        for (term_idx, i) in enumerate(indices)
            if checkbounds(Bool, dist_mut, matrix_idx, i) &&
               checkbounds(Bool, dist_wt, matrix_idx, i) &&
               checkbounds(Bool, PAE_mut, matrix_idx, i) &&
               checkbounds(Bool, PAE_wt, matrix_idx, i)
                
                d_mut = dist_mut[matrix_idx, i]
                d_wt = dist_wt[matrix_idx, i]

                if d_mut > 0.0 && d_wt > 0.0
                    mut_terms[round_idx, term_idx] = PAE_mut[matrix_idx, i]^(2 - rho) / (d_mut^α)
                    wt_terms[round_idx, term_idx] = PAE_wt[matrix_idx, i]^(2 - rho) / (d_wt^α)
                end
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
    α::Float64,
    config::EntropyConfig;
    verbose::Bool = false,
    cache_manager::Union{CacheManager,Nothing} = nothing,
)::Tuple{Vector{Float64},Vector{Float64},Vector{Float64}}
    coordinates = read_coordinates(pdb_path)
    Γ = compute_Γ(coordinates, config)
    mutations = read_mutations_from_file(task_file_path)
    results = Vector{Tuple{Float64,Float64,Float64}}()

    for m in mutations
        position = parse_mutation_position(m)
        result = process_single_mutation(m, position, single_ddG, paes, Γ, WT_pae, ddG_exp, rho, A, α, config, verbose; cache_manager)
        
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
    α::Float64,
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
    
    ΔΔS_val = ΔΔS(position, rho, α, Γ, pae_rounds, WT_pae, mutation, config; cache_manager)
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

"""
    optimize_gamma_parameters(coordinates, experimental_msf; verbose=false)

Optimizes Gamma kernel function parameters {η₁, κ₁, η₂, κ₂} to maximize Pearson correlation 
coefficient between theoretical and experimental MSF values.

# Arguments
- `coordinates::AbstractVector{<:AbstractVector{<:Real}}`: A vector of 3D coordinates for each atom.
- `experimental_msf::Vector{Float64}`: Experimental MSF (B-factor) values.
- `verbose::Bool`: Whether to print progress information.

# Returns
- `NamedTuple`: Best parameters and maximum PCC value.
"""
function optimize_gamma_parameters(
    coordinates::AbstractVector{<:AbstractVector{<:Real}},
    experimental_msf::Vector{Float64};
    verbose::Bool = false
)::NamedTuple{(:eta1, :kappa1, :eta2, :kappa2, :max_pcc), NTuple{5,Float64}}
    
    # Parameter ranges
    eta1_range = 1:30
    kappa1_range = 1:12
    eta2_range = 1:20
    kappa2_range = 1:12
    
    max_pcc = 0.0
    best_params = (eta1=1.0, kappa1=1.0, eta2=1.0, kappa2=1.0, max_pcc=0.0)
    
    total_combinations = length(eta1_range) * length(kappa1_range) * length(eta2_range) * length(kappa2_range)
    current_iteration = 0
    
    verbose && println("Starting parameter optimization with $total_combinations combinations...")
    
    for eta1 in eta1_range, kappa1 in kappa1_range, eta2 in eta2_range, kappa2 in kappa2_range
        current_iteration += 1
        
        # Calculate Gamma matrices
        Γ1 = Γ(coordinates, eta1, kappa1)
        Γ2 = Γ(coordinates, eta2, kappa2)
        Γ_total = Γ1 + Γ2
        
        # Calculate theoretical MSF
        theoretical_msf = msf(Γ_total)
        
        # Calculate Pearson correlation coefficient
        pcc = cor(theoretical_msf, experimental_msf)
        
        # Update best parameters if PCC improved
        if pcc > max_pcc
            max_pcc = pcc
            best_params = (eta1=Float64(eta1), kappa1=Float64(kappa1), 
                          eta2=Float64(eta2), kappa2=Float64(kappa2), max_pcc=max_pcc)
            verbose && println("New best PCC: $max_pcc with params η₁=$eta1, κ₁=$kappa1, η₂=$eta2, κ₂=$kappa2")
        end
        
        # Progress update
        if verbose && current_iteration % 1000 == 0
            progress = round(current_iteration / total_combinations * 100, digits=1)
            println("Progress: $progress% ($current_iteration/$total_combinations)")
        end
    end
    
    verbose && println("Optimization complete. Best PCC: $max_pcc")
    return best_params
end

"""
    fit_gamma_to_bfactor(pdb_file_path; verbose=false)

Convenience function to fit Gamma kernel parameters using B-factor data from a PDB file.

# Arguments
- `pdb_file_path::String`: Path to the PDB file.
- `verbose::Bool`: Whether to print progress information.

# Returns
- `NamedTuple`: Best parameters and maximum PCC value.
"""
function fit_gamma_to_bfactor(pdb_file_path::String; verbose::Bool = false)::NamedTuple
    # Read coordinates and B-factors
    coordinates = read_coordinates(pdb_file_path)
    b_factors_dict = get_residue_calpha_b_factor(pdb_file_path)
    
    verbose && println("Found $(length(b_factors_dict)) residues with B-factor data")
    
    # Extract B-factor values in residue order
    experimental_msf = Float64[]
    structure = read(pdb_file_path, PDBFormat)
    
    for residue in collectresidues(structure, standardselector)
        residue_id = "$(resname(residue)) $(resnumber(residue)) $(chainid(residue))"
        if haskey(b_factors_dict, residue_id)
            push!(experimental_msf, b_factors_dict[residue_id])
        else
            verbose && println("No B-factor found for residue: $residue_id")
        end
    end
    
    if length(experimental_msf) != length(coordinates)
        error("Mismatch between number of coordinates ($(length(coordinates))) and B-factors ($(length(experimental_msf))). Check if B-factor data is available in the PDB file.")
    end
    
    verbose && println("Loaded $(length(experimental_msf)) residues for optimization.")
    
    return optimize_gamma_parameters(coordinates, experimental_msf; verbose=verbose)
end

"""
    create_config_from_optimization(datadir, optimization_result; kwargs...)

Create an EntropyConfig from gamma parameter optimization results.

# Arguments
- `datadir::String`: Data directory path.
- `optimization_result::NamedTuple`: Result from optimize_gamma_parameters.
- `kwargs...`: Additional arguments for EntropyConfig constructor.

# Returns
- `EntropyConfig`: Configuration with optimized gamma parameters.
"""
function create_config_from_optimization(
    datadir::String,
    optimization_result::NamedTuple;
    kwargs...
)::EntropyConfig
    return EntropyConfig(
        datadir;
        eta1=round(Int, optimization_result.eta1),
        kappa1=round(Int, optimization_result.kappa1),
        eta2=round(Int, optimization_result.eta2),
        kappa2=round(Int, optimization_result.kappa2),
        kwargs...
    )
end
