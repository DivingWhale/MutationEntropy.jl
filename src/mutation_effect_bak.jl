"""
    ca_distMap(struc)

Calculate distance map between C-alpha atoms in the given protein structure.
"""
function ca_distMap(struc)
    dists = DistanceMap(collectatoms(struc, calphaselector))
    return dists
end

# Global cache for structure dimensions
const STRUCTURE_DIMENSION_CACHE = Dict{String, Int}()

"""
    get_structure_dimension(datadir::String, mutation::String, round::Int; use_cache::Bool=true)

Get the dimension of protein structure (number of C-alpha atoms) for initializing distance matrix.
Supports caching to avoid repeated calculations for the same protein variant.

# Arguments
- `datadir::String`: Base directory containing protein data
- `mutation::String`: Mutation identifier (e.g., "A140E")
- `round::Int`: Round number for the simulation
- `use_cache::Bool=true`: Whether to use cached dimension values

# Returns
- `Int`: Number of C-alpha atoms in the structure
"""
function get_structure_dimension(datadir::String, mutation::String, round::Int; use_cache::Bool=true)
    # Create a unique key for this protein variant
    cache_key = mutation
    
    # Check if we have a cached value
    if use_cache && haskey(STRUCTURE_DIMENSION_CACHE, cache_key)
        return STRUCTURE_DIMENSION_CACHE[cache_key]
    end
    
    # If not in cache, calculate the dimension
    subfolder_pattern = "$(datadir)/$(mutation)_$(round)/seed-*_sample-0"
    subfolders = glob(subfolder_pattern)
    
    if isempty(subfolders)
        error("No matching folder found: $subfolder_pattern")
    end
    
    subfolder = subfolders[1]
    model_pdb = joinpath(subfolder, "model.pdb")
    
    if !isfile(model_pdb)
        # Convert CIF to PDB using Python script
        cmd = `python src/convert.py $subfolder`
        run(cmd)
        
        if !isfile(model_pdb)
            error("Unable to create PDB file: $model_pdb")
        end
    end
    
    # Read structure and get dimension
    structure = read(model_pdb, PDBFormat)
    ca_atoms = collectatoms(structure, calphaselector)
    dim = length(ca_atoms)
    
    # Cache the value for future use
    if use_cache
        STRUCTURE_DIMENSION_CACHE[cache_key] = dim
    end
    
    return dim
end

"""
    clear_dimension_cache()

Clear the cached protein structure dimensions.
"""
function clear_dimension_cache()
    empty!(STRUCTURE_DIMENSION_CACHE)
    return nothing
end

"""
    get_dimension_cache()

Get a copy of the current dimension cache.
"""
function get_dimension_cache()
    return copy(STRUCTURE_DIMENSION_CACHE)
end

"""
    get_dist_map(datadir::String, mutation::String, round::Int; use_cache::Bool=true)

Calculate the distance map for a specific mutation and round.
Returns a matrix where each element represents the distance between C-alpha atoms.

# Arguments
- `datadir::String`: Base directory containing protein data
- `mutation::String`: Mutation identifier (e.g., "A140E")
- `round::Int`: Round number to process
- `use_cache::Bool=true`: Whether to use cached dimension values

# Returns
- `Matrix{Float64}`: Distance matrix
"""
function get_dist_map(datadir::String, mutation::String, round::Int; use_cache::Bool=true)
    
    subfolder_pattern = "$(datadir)/$(mutation)_$(round)/seed-*_sample-0"
    subfolders = glob(subfolder_pattern)
    
    if isempty(subfolders)
        error("No matching folder found: $subfolder_pattern")
    end
    
    subfolder = subfolders[1]
    model_pdb = joinpath(subfolder, "model.pdb")
    
    if !isfile(model_pdb)
        # Convert CIF to PDB using Python script
        cmd = `python src/convert.py $subfolder`
        run(cmd)
        
        if !isfile(model_pdb)
            error("Unable to create PDB file: $model_pdb")
        end
    end
    
    d_matrix_file = joinpath(subfolder, "d_matrix.jld2")
    
    if !isfile(d_matrix_file)
        # Calculate the distance matrix
        structure = read(model_pdb, PDBFormat)
        d_matrix = ca_distMap(structure).data
        
        # Save the distance matrix using Julia's serialization
        jldsave(d_matrix_file; d_matrix=d_matrix)
    else
        # Load the precomputed distance matrix
        d_matrix = load(d_matrix_file, "d_matrix")
    end
    
    return d_matrix
end

"""
    precompute_dimensions(datadir::String, mutations::Vector{String}, round::Int=1)

Precompute and cache dimensions for a list of mutations to speed up subsequent processing.
This is useful when processing multiple mutations of the same protein.

# Arguments
- `datadir::String`: Base directory containing protein data
- `mutations::Vector{String}`: List of mutation identifiers (e.g., ["A140E", "A140G"])
- `round::Int=1`: The round number to use for dimension calculation

# Returns
- `Dict{String, Int}`: A dictionary mapping mutations to their dimensions
"""
function precompute_dimensions(datadir::String, mutations::Vector{String}, round::Int=1)
    dimensions = Dict{String, Int}()
    
    for mutation in mutations
        try
            dim = get_structure_dimension(datadir, mutation, round, use_cache=true)
            dimensions[mutation] = dim
            println("Cached dimension for $mutation: $dim")
        catch e
            @warn "Failed to compute dimension for $mutation: $e"
        end
    end
    
    return dimensions
end

"""
    process_multiple_mutations(datadir::String, mutations::Vector{String}, round::Int; 
                               precache::Bool=true, parallel::Bool=false)

Process multiple mutations and calculate distance maps for each for a specific round.
Supports precomputing dimensions and optional parallel processing.

# Arguments
- `datadir::String`: Base directory containing protein data
- `mutations::Vector{String}`: List of mutation identifiers (e.g., ["A140E", "A140G"])
- `round::Int`: Round number to process
- `precache::Bool=true`: Whether to precache dimensions before processing
- `parallel::Bool=false`: Whether to use parallel processing (requires Distributed package)

# Returns
- `Dict{String, Matrix{Float64}}`: A dictionary mapping mutations to their distance matrices
"""
function process_multiple_mutations(datadir::String, mutations::Vector{String}, 
                                   round::Int; 
                                   precache::Bool=true, 
                                   parallel::Bool=false)
    
    # Optional precomputing of dimensions
    if precache
        precompute_dimensions(datadir, mutations, round)
    end
    
    results = Dict{String, Matrix{Float64}}()
    
    if parallel
        # Import Distributed module functions if parallel is enabled
        # This ensures remotecall_fetch is available
        results = Dict(
            mutation => remotecall_fetch(
                get_dist_map, 
                workers()[mod1(i, nworkers())], 
                datadir, mutation, round
            ) for (i, mutation) in enumerate(mutations)
        )
    else
        # Sequential processing
        for mutation in mutations
            results[mutation] = get_dist_map(datadir, mutation, round)
            println("Processed $mutation")
        end
    end
    
    return results
end

"""
    save_dimension_cache(filename::String)

Save the current dimension cache to a file.

# Arguments
- `filename::String`: Path to the file where cache will be saved
"""
function save_dimension_cache(filename::String)
    jldsave(filename; cache=STRUCTURE_DIMENSION_CACHE)
    return nothing
end

"""
    load_dimension_cache(filename::String)

Load dimension cache from a file.

# Arguments
- `filename::String`: Path to the file containing saved cache
"""
function load_dimension_cache(filename::String)
    if isfile(filename)
        data = load(filename)
        merge!(STRUCTURE_DIMENSION_CACHE, data["cache"])
    else
        error("Cache file not found: $filename")
    end
    return nothing
end

"""
    get_low_plddt_residues(mutation::String, round::Int, datadir::String; threshold::Float64=90.0)

Identify residues with low pLDDT scores (predicted Local Distance Difference Test) from AlphaFold models.
Lower pLDDT scores indicate lower confidence regions in the predicted structure.

# Arguments
- `mutation::String`: Mutation identifier (e.g., "A140E")
- `round::Int`: Round number to analyze
- `datadir::String`: Base directory containing protein data
- `threshold::Float64=90.0`: pLDDT threshold below which residues are considered low confidence

# Returns
- `Vector{Int}`: List of residue IDs with pLDDT scores below the threshold
"""
function get_low_plddt_residues(mutation::String, round::Int, datadir::String; threshold::Float64=90.0)
    subfolder_pattern = "$(datadir)/$(mutation)_$(round)/seed-*_sample-0"
    subfolders = glob(subfolder_pattern)
    
    if isempty(subfolders)
        error("No matching folder found: $subfolder_pattern")
    end
    
    subfolder = subfolders[1]
    model_pdb = joinpath(subfolder, "model.pdb")
    temp_pdb = "$(mutation)_temp_pl.pdb"
    
    # Check if model.pdb exists, if not create it from model.cif
    if !isfile(model_pdb)
        cmd = `python src/convert.py $subfolder`
        run(cmd)
        
        if !isfile(model_pdb)
            error("Unable to create PDB file: $model_pdb")
        end
    end
    
    # Load structure using BioStructures
    structure = read(model_pdb, PDBFormat)
    
    # Load confidence scores from JSON file
    confidence_file = joinpath(subfolder, "confidences.json")
    if !isfile(confidence_file)
        error("Confidence file not found: $confidence_file")
    end
    
    confidence_data = JSON.parsefile(confidence_file)
    atom_plddts = confidence_data["atom_plddts"]
    
    # Process residue-level pLDDT scores
    residue_plddts = Dict{Int, Float64}()
    atom_index = 1
    
    for model in structure
        for chain in model
            for residue in chain
                res_id = BioStructures.resid(residue) |> x -> parse(Int, x)
                residue_atoms = collect(atoms(residue))
                
                if !haskey(residue_plddts, res_id) && !isempty(residue_atoms)
                    # Calculate mean pLDDT for atoms in this residue
                    atom_count = length(residue_atoms)
                    
                    # Ensure we have enough pLDDT values
                    if atom_index + atom_count - 1 <= length(atom_plddts)
                        residue_plddt = mean(atom_plddts[atom_index:(atom_index + atom_count - 1)])
                        residue_plddts[res_id] = residue_plddt
                    end
                    
                    atom_index += atom_count
                end
            end
        end
    end
    
    # Identify residues below threshold
    low_plddt_residues = [res_id for (res_id, plddt) in residue_plddts if plddt < threshold]
    
    println("Found $(length(low_plddt_residues)) residues with pLDDT < $threshold for $mutation round $round")
    return sort(low_plddt_residues)
end

"""
    find_residues_within_distance(residue_id::Int, datadir::String, round::Int, 
                                 mutation::String; distance::Float64=13.0)

Find all residues within a specific distance of a target residue for a given round.

# Arguments
- `residue_id::Int`: The ID of the target residue
- `datadir::String`: Base directory containing protein data
- `round::Int`: Round number to process
- `mutation::String`: Mutation identifier (e.g., "A140E")
- `distance::Float64=13.0`: Distance cutoff in Angstroms

# Returns
- `Vector{Int}`: List of residue IDs within the distance threshold
"""
function find_residues_within_distance(residue_id::Int, d_matrix::Matrix{Float64}; distance::Float64=13.0)
    
    # Find residues within the distance threshold
    nearby_residues = findall(x -> x < distance, d_matrix[residue_id, :])
    
    # println("Found $(length(nearby_residues)) residues within $(distance)Å of residue $residue_id")
    
    return sort(nearby_residues)
end

"""
    read_pae(datadir::String, mutation::String, round::Int; use_cache::Bool=true)

Read the Predicted Aligned Error (PAE) matrix for a specific round of AlphaFold prediction.
PAE values represent AlphaFold's estimate of position error between residue pairs.
Uses dimension caching to validate matrix dimensions.

# Arguments
- `datadir::String`: Base directory containing protein data
- `mutation::String`: Mutation identifier (e.g., "A140E")
- `round::Int`: Round number to process
- `use_cache::Bool=true`: Whether to use the cached dimension values

# Returns
- `Matrix{Float64}`: PAE matrix for the specified round
"""
function read_paes(datadir::String, mutation::String, round::Int; use_cache::Bool=true)
    # Find the appropriate subfolder
    subfolder_pattern = "$(datadir)/$(mutation)_$(round)/seed-*_sample-0"
    subfolders = glob(subfolder_pattern)
    
    if isempty(subfolders)
        error("No matching folder found for $mutation round $round")
    end
    
    subfolder = subfolders[1]
    confidence_file = joinpath(subfolder, "confidences.json")
    
    if !isfile(confidence_file)
        error("Confidence file not found for $mutation round $round: $confidence_file")
    end
    
    # Load PAE matrix from the confidence file
    try
        # Read and parse the PAE matrix
        raw_data = JSON.parsefile(confidence_file)
        pae_any = raw_data["pae"]
        
        # Convert all elements to Float64 (similar to read_single_pae)
        pae = map(x -> map(Float64, x), pae_any)
        
        # Convert to matrix format and ensure it's a standard Matrix{Float64}
        pae = Matrix(hcat(pae...)')
        
        # Verify dimensions if using cache
        if use_cache
            dim = get_structure_dimension(datadir, mutation, round, use_cache=true)
            if size(pae) != (dim, dim)
                @warn "PAE matrix dimension mismatch for $mutation round $round: expected ($dim, $dim), got $(size(pae))."
            end
        end
        
        println("Processed PAE matrix for $mutation round $round")
        return pae
        
    catch e
        error("Error processing PAE data for $mutation round $round: $e")
    end
end

function calculate_strain(alpha::Float64, paes::Matrix{Float64}, d_matrix::Matrix{Float64}, mutation::String)
    r_index = parse(Int, match(r"[A-Z](\d+)[A-Z]", mutation).captures[1])
    nearby_residues = find_residues_within_distance(r_index, d_matrix)
    strain = 0.0
    for residue in nearby_residues
        strain += paes[r_index, residue] / (d_matrix[r_index, residue]^alpha)
    end
    strain = strain / length(nearby_residues)

    return strain
end