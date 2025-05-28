# filepath: /hpc2hdd/home/sxu128/MD/MutationEntropy/src/mutation_effect.jl
"""
    ca_distMap(struc)

Calculate distance map between C-alpha atoms in the given protein structure.
"""
function ca_distMap(struc)
    dists = DistanceMap(collectatoms(struc, calphaselector))
    return dists
end

"""
    get_structure_dimension(datadir::String, mutation::String, round::Int)

Get the dimension of protein structure (number of C-alpha atoms) for initializing distance matrix.

# Arguments
- `datadir::String`: Base directory containing protein data
- `mutation::String`: Mutation identifier (e.g., "A140E")
- `round::Int`: Round number for the simulation

# Returns
- `Int`: Number of C-alpha atoms in the structure
"""
function get_structure_dimension(datadir::String, mutation::String, round::Int)
    # Construct the path to the model PDB file
    subfolder_pattern = "$(datadir)/$(mutation)_$(round)/seed-*_sample-0"
    subfolders = glob(subfolder_pattern)
    
    if isempty(subfolders)
        error("No matching folder found: $subfolder_pattern")
    end
    
    subfolder = subfolders[1]
    model_pdb = joinpath(subfolder, "model.pdb")
    
    if !isfile(model_pdb)
        # Convert CIF to PDB using Python script with conda environment
        script_path = abspath(joinpath(@__DIR__, "convert.py"))
        subfolder_abs = abspath(subfolder)
        conda_env = "bio"
        
        cmd = `bash -c "source ~/.bashrc && conda activate $conda_env && python $script_path $subfolder_abs"`
        println("Running conversion with conda: $cmd")
        run(cmd)
        
        if !isfile(model_pdb)
            error("Unable to create PDB file: $model_pdb")
        end
    end
    
    # Read structure and get dimension
    structure = read(model_pdb, PDBFormat)
    ca_atoms = collectatoms(structure, calphaselector)
    dim = length(ca_atoms)
    
    return dim
end

"""
    get_dist_map(datadir::String, mutation::String, round::Int)

Calculate the distance map for a specific mutation and round.
Returns a matrix where each element represents the distance between C-alpha atoms.

# Arguments
- `datadir::String`: Base directory containing protein data
- `mutation::String`: Mutation identifier (e.g., "A140E")
- `round::Int`: Round number to process

# Returns
- `Matrix{Float64}`: Distance matrix
"""
function get_dist_map(datadir::String, mutation::String, round::Int)
    
    subfolder_pattern = "$(datadir)/$(mutation)_$(round)/seed-*_sample-0"
    subfolders = glob(subfolder_pattern)
    
    if isempty(subfolders)
        error("No matching folder found: $subfolder_pattern")
    end
    
    subfolder = subfolders[1]
    model_pdb = joinpath(subfolder, "model.pdb")
    
    if !isfile(model_pdb)
        # Check if CIF file exists
        cif_file = joinpath(subfolder, "model.cif")
        if !isfile(cif_file)
            error("Neither PDB nor CIF file found in $subfolder")
        end

        # Get absolute paths to ensure correct file access
        script_path = abspath(joinpath(@__DIR__, "convert.py"))
        subfolder_abs = abspath(subfolder)
        
        try
            # Use the conda environment by activating it first
            conda_env = "bio"  # Your conda environment name
            
            # Create a command that activates the conda environment and runs the Python script
            cmd = `bash -c "source ~/.bashrc && conda activate $conda_env && python $script_path $subfolder_abs"`
            run(cmd)
            
            if !isfile(model_pdb)
                error("Conversion completed but PDB file not created at $model_pdb")
            end
        catch e
            println("Error during conversion: $e")
            error("Failed to convert CIF to PDB. See above error message.")
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

Precompute dimensions for a list of mutations to speed up subsequent processing.
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
            dim = get_structure_dimension(datadir, mutation, round)
            dimensions[mutation] = dim
            println("Calculated dimension for $mutation: $dim")
        catch e
            @warn "Failed to compute dimension for $mutation: $e"
        end
    end
    
    return dimensions
end

"""
    process_multiple_mutations(datadir::String, mutations::Vector{String}, round::Int; 
                               parallel::Bool=false)

Process multiple mutations and calculate distance maps for each for a specific round.
Supports optional parallel processing.

# Arguments
- `datadir::String`: Base directory containing protein data
- `mutations::Vector{String}`: List of mutation identifiers (e.g., ["A140E", "A140G"])
- `round::Int`: Round number to process
- `parallel::Bool=false`: Whether to use parallel processing (requires Distributed package)

# Returns
- `Dict{String, Matrix{Float64}}`: A dictionary mapping mutations to their distance matrices
"""
function process_multiple_mutations(datadir::String, mutations::Vector{String}, 
                                   round::Int;
                                   parallel::Bool=false)
    
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
        # Convert CIF to PDB using Python script with conda environment
        script_path = abspath(joinpath(@__DIR__, "convert.py"))
        subfolder_abs = abspath(subfolder)
        conda_env = "bio"
        
        cmd = `bash -c "source ~/.bashrc && conda activate $conda_env && python $script_path $subfolder_abs"`
        println("Running conversion with conda: $cmd")
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
    
    return sort(low_plddt_residues)
end

"""
    find_residues_within_distance(residue_id::Int, d_matrix::Matrix{Float64}; distance::Float64=13.0)

Find all residues within a specific distance of a target residue.

# Arguments
- `residue_id::Int`: The ID of the target residue
- `d_matrix::Matrix{Float64}`: Distance matrix containing pairwise distances between residues
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
    read_pae(datadir::String, mutation::String, round::Int)

Read the Predicted Aligned Error (PAE) matrix for a specific round of AlphaFold prediction.
PAE values represent AlphaFold's estimate of position error between residue pairs.

# Arguments
- `datadir::String`: Base directory containing protein data
- `mutation::String`: Mutation identifier (e.g., "A140E")
- `round::Int`: Round number to process

# Returns
- `Matrix{Float64}`: PAE matrix for the specified round
"""
function read_pae(datadir::String, mutation::String, round::Int)
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
        
        # Convert all elements to Float64
        pae = map(x -> map(Float64, x), pae_any)
        
        # Convert to matrix format
        pae = Matrix(hcat(pae...)')
        
        return pae
        
    catch e
        error("Error processing PAE data for $mutation round $round: $e")
    end
end

function calculate_strain(alpha::Float64, pae::Matrix{Float64}, d_matrix::Matrix{Float64}, site::Int64)
    nearby_residues = find_residues_within_distance(site, d_matrix)
    strain = 0.0
    for residue in nearby_residues
        if residue == site
            continue
        end
        strain += pae[site, residue] / (d_matrix[site, residue]^alpha)
    end
    strain = strain / length(nearby_residues)

    return strain
end

"""
    collect_strains(datadir::String, protein::String; 
                    alpha::Float64=2.0, 
                    residue_range=nothing, 
                    num_rounds::Int=20, 
                    verbose::Bool=true,
                    cache::Bool=true,
                    plddt_threshold::Float64=90.0)

Collect and average strain values for multiple residues across multiple rounds.
Low pLDDT residues are automatically excluded from calculations by default.

# Arguments
- `datadir::String`: Base directory containing protein data
- `protein::String`: Protein identifier (e.g., "thermonulease")
- `alpha::Float64=2.0`: Exponent parameter for strain calculation
- `residue_range=nothing`: Range of residues to analyze (e.g., 83:231), if nothing analyzes all residues
- `num_rounds::Int=20`: Number of rounds to average over
- `verbose::Bool=true`: Whether to print progress information
- `cache::Bool=true`: Whether to cache distance maps and PAE matrices
- `plddt_threshold::Float64=90.0`: pLDDT threshold below which residues are considered low confidence

# Returns
- `Dict{Int, Float64}`: Dictionary mapping residue IDs to their averaged strain values
"""
function collect_strains(datadir::String, protein::String; 
                         alpha::Float64=2.0, 
                         residue_range=nothing, 
                         num_rounds::Int=20,
                         verbose::Bool=true,
                         cache::Bool=true,
                         plddt_threshold::Float64=90.0)
    
    # Determine residue range if not provided
    if residue_range === nothing
        # Use round 1 to get the dimension of the protein
        dists_r1 = get_dist_map(datadir, protein, 1)
        residue_range = 1:size(dists_r1, 1)
        if verbose
            println("Residue range not provided, using full range: $(first(residue_range)) to $(last(residue_range))")
        end
    end
    
    # Filter out low pLDDT residues (now default behavior)
    try
        low_plddt_residues = get_low_plddt_residues(protein, 1, datadir, threshold=plddt_threshold)
        if !isempty(low_plddt_residues)
            # Create a new range excluding low pLDDT residues
            filtered_residue_range = filter(r -> !(r in low_plddt_residues), residue_range)
            
            if verbose
                excluded_count = length(residue_range) - length(filtered_residue_range)
                println("Excluded $excluded_count residues with pLDDT < $plddt_threshold")
            end
            
            residue_range = filtered_residue_range
        end
    catch e
        @warn "Failed to filter low pLDDT residues: $e"
    end
    
    # Create cache directory if needed
    cache_dir = joinpath(datadir, "strain_cache")
    if cache && !isdir(cache_dir)
        mkpath(cache_dir)
    end
    
    # Initialize a dictionary with keys from the residue range
    site_value_sums = Dict(k => 0.0 for k in residue_range)
    
    # Pre-compute nearby residues for each site (this avoids redundant calculations)
    nearby_residues_cache = Dict{Int, Vector{Int}}()
    
    # Process rounds
    function process_round(r)
        local_results = Dict(k => 0.0 for k in residue_range)
        
        # Check cache first
        # Include information about pLDDT filtering in the cache file name
        plddt_suffix = "_plddt$(plddt_threshold)"
        cache_file = joinpath(cache_dir, "$(protein)_round$(r)_alpha$(alpha)$(plddt_suffix).jld2")
        
        if cache && isfile(cache_file)
            if verbose
                println("Loading cached data for round $r")
            end
            try
                local_results = load(cache_file, "strain_data")
                return local_results
            catch e
                @warn "Failed to load cached data for round $r: $e. Recomputing..."
            end
        end
        
        try
            # Get distance map and PAE matrix for this round
            dists = get_dist_map(datadir, protein, r)
            paes = read_pae(datadir, protein, r)
            
            # Calculate strain for each site in the range
            for site in residue_range
                # Get or compute nearby residues
                if !haskey(nearby_residues_cache, site)
                    nearby_residues_cache[site] = find_residues_within_distance(site, dists)
                end
                
                nearby = nearby_residues_cache[site]
                strain_sum = 0.0
                valid_count = 0
                
                # Calculate strain for this site
                for residue in nearby
                    if residue == site
                        continue
                    end
                    
                    # Ensure we're within bounds
                    if site <= size(paes, 1) && residue <= size(paes, 2) && 
                       site <= size(dists, 1) && residue <= size(dists, 2)
                        strain_sum += paes[site, residue] / (dists[site, residue]^alpha)
                        valid_count += 1
                    end
                end
                
                # Avoid division by zero
                if valid_count > 0
                    local_results[site] = strain_sum / valid_count
                else
                    local_results[site] = 0.0
                    @warn "No valid neighbors found for site $site in round $r"
                end
            end
            
            # Cache the results
            if cache
                try
                    jldsave(cache_file; strain_data=local_results)
                catch e
                    @warn "Failed to cache results for round $r: $e"
                end
            end
            
            if verbose
                println("Completed round $r")
            end
            
            return local_results
            
        catch e
            @error "Failed to process round $r: $e"
            return Dict(k => 0.0 for k in residue_range)  # Return zeros on error
        end
    end
    
    # Process rounds sequentially
    for r in 1:num_rounds
        if verbose
            println("Processing round $r of $num_rounds...")
        end
        
        r_results = process_round(r)
        
        # Aggregate results
        for (site, value) in r_results
            site_value_sums[site] += value
        end
    end
    
    # Calculate the average values
    final_averaged_values = Dict{eltype(residue_range), Float64}()
    for site in residue_range
        if haskey(site_value_sums, site)
            final_averaged_values[site] = site_value_sums[site] / num_rounds
        else
            @warn "Site $site was not found in site_value_sums. This should not happen."
        end
    end
    
    # No save functionality
    
    # Report results
    if verbose
        println("\nProcessing completed")
        println("\nFinal Averaged Values per Site:")
        for (site, avg_value) in sort(collect(final_averaged_values))
            println("Site $site: $avg_value")
        end
    end
    
    return final_averaged_values
end

function calculate_ME(Strain::Dict{Int, Float64}, Strain_wt::Dict{Int, Float64})
    # Find the intersection of keys - residues present in both dictionaries
    common_sites = intersect(keys(Strain), keys(Strain_wt))
    
    if isempty(common_sites)
        error("No common residue sites found between the two strain dictionaries")
    end
    
    # Calculate the Mutation Entropy (ME) as the sum of differences for common sites
    ME = Dict{Int, Float64}()
    for site in common_sites
        ME[site] = abs(Strain[site] - Strain_wt[site])
    end
    
    # Log information about the intersection
    missing_count = length(keys(Strain)) + length(keys(Strain_wt)) - 2*length(common_sites)
    if missing_count > 0
        @info "$(missing_count) residue sites were excluded from ME calculation due to missing in either strain dataset"
        @info "Calculated ME for $(length(common_sites)) common residue sites"
    end
    
    return ME
end