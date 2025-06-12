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

function read_single_pae(data_path::String, mutation::String)
    pae_files = glob("$(mutation)_*", data_path)
    pae_data = nothing
    num_files = 0
    for pae_file in pae_files
        top1_files = glob("seed-*_sample-0", pae_file)
        if !isempty(top1_files)
            top1 = top1_files[1]
            raw_data = JSON.parsefile(joinpath(top1, "confidences.json"))
            pae_any = raw_data["pae"]
            pae = map(x -> map(Float64, x), pae_any)
            pae =  hcat(pae...)'
          
            if pae_data === nothing
                pae_data = zeros(Float64, size(pae))
            end
           
            pae_data .+= pae
            num_files += 1
        end
    end
    if num_files > 0
        return pae_data ./ num_files
    else
        return nothing
    end
end

function read_paes(task_file_path::String, source_data_path::String)
    CACHE_DIR = ".pae_cache" # Cache directory will be created in the current working directory
    println("Initializing PAE reading. Cache directory: $(abspath(CACHE_DIR))")

    try
        mkpath(CACHE_DIR)
        println("Cache directory ensured at $(abspath(CACHE_DIR))")
    catch e
        @warn "Could not create cache directory $CACHE_DIR: $e. Caching may not work as expected."
    end

    mutations_list = read_mutations_from_file(task_file_path)
    num_mutations = length(mutations_list)
    println("Found $(num_mutations) mutations to process from $task_file_path.")
    
    results_channel = Channel{Tuple{String, Union{Matrix{Float64}, Nothing}}}(num_mutations)

    println("Starting parallel processing of mutations...")
    processed_count_atomic = Threads.Atomic{Int}(0) # Atomic counter for progress

    Threads.@threads for mutation_original_case in mutations_list
        mutation = string(strip(mutation_original_case))
        
        current_processed_count = Threads.atomic_add!(processed_count_atomic, 1)
        println("Processing $(current_processed_count)/$(num_mutations): $mutation")

        cache_file_name = "$(mutation).jld2"
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
            loaded_pae = read_single_pae(source_data_path, mutation)
            
            if loaded_pae !== nothing
                try
                    JLD2.save_object(cache_file_path, loaded_pae)
                catch err
                    @warn "Failed to save $mutation to cache file $cache_file_path: $err."
                end
            else
                println("No PAE data found for $mutation from source $source_data_path.") # Kept as it's an important specific outcome
            end
        end
        
        put!(results_channel, (mutation, loaded_pae))
    end
    
    close(results_channel)
    println("Finished parallel processing. Aggregating results...")

    paes = Dict{String, Matrix{Float64}}()
    aggregated_count = 0
    for (mut, pae_matrix) in results_channel
        if pae_matrix !== nothing
            paes[mut] = pae_matrix
        else
            println("Warning: No PAE data for mutation $mut after processing.")
        end
        aggregated_count += 1
    end
    # The variable `processed_count` was used in the previous version for the final print.
    # Using `aggregated_count` which is the actual number of items taken from the channel.
    # Or, using `num_mutations` to reflect total attempted.
    # Let's use `num_mutations` for total attempted and `length(paes)` for successful.
    println("Finished aggregation. Processed $(num_mutations) attempted mutations. $(length(paes)) mutations have PAE data.")
    
    return paes
end

function ΔΔS(position::Int64, rho::Float64, Γ::Matrix{Float64}, PAE_mut::Matrix{Float64}, PAE_wt::Matrix{Float64}, offset::Int64=0)
    matrix_idx = position - offset
    indices = findall(x -> isapprox(x, 0.0, atol=1e-6), Γ[matrix_idx, :])
    ΔΔS = 0.0
    for i in indices
        ΔΔS += PAE_mut[matrix_idx, i]^(2-rho) - PAE_wt[matrix_idx, i]^(2-rho)
    end
    ΔΔS = ΔΔS ./ length(indices)
    return ΔΔS
end

function ΔΔG_prime(A::Float64, ΔΔS::Float64, ΔΔG::Float64)
    return ΔΔG + A * ΔΔS
end

"""
    calculate_ddgs(task_file_path::String, single_ddG::Dict{String, Float64}, pdb_path::String, WT_pae::Matrix{Float64}, paes::Dict{String, Matrix{Float64}}, ddG_exp::DataFrame, rho::Float64, A::Float64, offset::Int64=0, verbose::Bool=false)

Calculate the predicted ΔΔG values for a set of mutations.

# Arguments
- `task_file_path::String`: Path to file containing list of mutations to analyze
- `single_ddG::Dict{String, Float64}`: Dictionary mapping mutation strings to Rosetta ddG values
- `pdb_path::String`: Path to the PDB file for computing the Gamma matrix
- `WT_pae::Matrix{Float64}`: PAE matrix for wild-type protein
- `paes::Dict{String, Matrix{Float64}}`: Dictionary mapping mutations to their PAE matrices
- `ddG_exp::DataFrame`: DataFrame containing experimental ddG values
- `rho::Float64`: Parameter controlling contribution of PAE differences
- `A::Float64`: Scaling parameter for entropy contribution
- `offset::Int64`: Offset between protein position numbers and matrix indices (default: 0)
- `verbose::Bool`: If true, prints details about skipped variants

# Returns
- `filtered_ddG_exp::Vector{Float64}`: Filtered experimental ddG values
- `ΔΔGs::Vector{Float64}`: Predicted total ddG values including entropy contribution
- `r_ddGs::Vector{Float64}`: Original Rosetta ddG predictions
"""
function calculate_ddgs(task_file_path::String, single_ddG::Dict{String, Float64}, pdb_path::String, WT_pae::Matrix{Float64}, paes::Dict{String, Matrix{Float64}}, ddG_exp::DataFrame, rho::Float64, A::Float64, offset::Int64=0, verbose::Bool=false)
    # Compute Gamma matrix from PDB file
    coordinates = read_coordinates(pdb_path)
    Γ = compute_Γ(coordinates)
    
    # Process mutations
    mutations = read_mutations_from_file(task_file_path)
    results = Vector{Tuple{Float64, Float64, Float64}}()
    
    for m in mutations
        position = parse_mutation_position(m)
        result = process_single_mutation(m, position, single_ddG, paes, Γ, WT_pae, ddG_exp, rho, A, offset, verbose)
        if result !== nothing && result[2] !== NaN
            push!(results, result)
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
                                paes::Dict{String, Matrix{Float64}}, Γ::Matrix{Float64}, 
                                WT_pae::Matrix{Float64}, ddG_exp::DataFrame, rho::Float64, A::Float64, offset::Int64=0, verbose::Bool=false)
    
    mutation_upper = uppercase(mutation)
    if !haskey(single_ddG, mutation_upper)
        if verbose
            println("Skipping variant: $mutation (not found in single_ddG)")
        end
        return nothing
    end
    
    ddG = single_ddG[mutation_upper]
    pae = paes[mutation]
    ΔΔS = MutationEntropy.ΔΔS(position, rho, Γ, pae, WT_pae, offset)
    predicted_ddG = ΔΔG_prime(A, ΔΔS, ddG)
    experimental_ddG = get_experimental_ddg(ddG_exp, position, string(mutation_upper[end]))
    
    return (experimental_ddG, predicted_ddG, ddG)
end