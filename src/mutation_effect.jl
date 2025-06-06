"""
    ca_distMap(struc)

Calculate distance map between C-alpha atoms in the given protein structure.
"""
function ca_distMap(struc) # Assuming struc is a ProteinStructure object
    # This function is already concise.
    # Ensure calphaselector is correctly defined/imported (likely from BioStructures).
    dists = DistanceMap(collectatoms(struc, calphaselector))
    return dists
end

"""
    _ensure_pdb_from_cif_if_needed(subfolder::String, model_pdb_path::String)

Helper to convert CIF to PDB using a Python script if model.pdb doesn't exist.
This is an internal helper to reduce code duplication in get_structure_dimension and get_dist_map.
"""
function _ensure_pdb_from_cif_if_needed(subfolder::String, model_pdb_path::String)
    if isfile(model_pdb_path)
        return true # PDB already exists
    end

    # Check if CIF file exists as a source for conversion
    cif_file = joinpath(subfolder, "model.cif")
    if !isfile(cif_file)
        # Error message tailored to where it might be called from get_dist_map
        # For get_structure_dimension, the error is slightly different in original.
        # This helper standardizes it. If a more specific error is needed,
        # this part might need to be duplicated or passed in.
        error("Neither PDB nor CIF file found in $subfolder") 
    end

    # Convert CIF to PDB using Python script with conda environment
    script_path = abspath(joinpath(@__DIR__, "convert.py"))
    subfolder_abs = abspath(subfolder)
    conda_env = "bio" # As used in the original code

    cmd = `bash -c "source ~/.bashrc && conda activate $conda_env && python $script_path $subfolder_abs"`
    
    try
        run(cmd)
    catch e
        # This combines error messages from original get_structure_dimension and get_dist_map
        println("Error during conversion script execution: $e") 
        error("Failed to convert CIF to PDB. Check conda environment '$conda_env', Python script '$script_path', and related files.")
    end
    
    if !isfile(model_pdb_path)
        # Error message consistency
        error("Conversion attempted but PDB file not created: $model_pdb_path")
    end
    return true
end


"""
    get_structure_dimension(datadir::String, mutation::String, round::Int)

Get the dimension of protein structure (number of C-alpha atoms) for initializing distance matrix.
"""
function get_structure_dimension(datadir::String, mutation::String, round_val::Int) # Renamed round to round_val to avoid conflict if `round` function is meant
    # Construct the path to the model PDB file
    subfolder_pattern = "$(datadir)/$(mutation)_$(round_val)/seed-*_sample-0"
    subfolders = glob(subfolder_pattern)
    
    if isempty(subfolders)
        error("No matching folder found: $subfolder_pattern")
    end
    
    subfolder = subfolders[1] # Assumes the first match is the desired one
    model_pdb = joinpath(subfolder, "model.pdb")
    
    _ensure_pdb_from_cif_if_needed(subfolder, model_pdb) # Handles PDB creation
    
    # Read structure and get dimension
    structure = read(model_pdb, PDBFormat)
    ca_atoms = collectatoms(structure, calphaselector)
    dim = length(ca_atoms)
    
    return dim
end

"""
    get_dist_map(datadir::String, mutation::String, round::Int)

Calculate the distance map for a specific mutation and round.
"""
function get_dist_map(datadir::String, mutation::String, round_val::Int) # Renamed round
    subfolder_pattern = "$(datadir)/$(mutation)_$(round_val)/seed-*_sample-0"
    subfolders = glob(subfolder_pattern)
    
    if isempty(subfolders)
        error("No matching folder found: $subfolder_pattern")
    end
    
    subfolder = subfolders[1]
    model_pdb = joinpath(subfolder, "model.pdb")
    
    _ensure_pdb_from_cif_if_needed(subfolder, model_pdb) # Handles PDB creation

    d_matrix_file = joinpath(subfolder, "d_matrix.jld2")
    
    local d_matrix # Ensure d_matrix is in scope
    if !isfile(d_matrix_file)
        structure = read(model_pdb, PDBFormat)
        d_matrix = ca_distMap(structure).data # .data extracts the matrix
        
        try
            jldsave(d_matrix_file; d_matrix=d_matrix)
        catch e
            @warn "Failed to save distance matrix to $d_matrix_file: $e" # Added warning for save failure
        end
    else
        try
            d_matrix = load(d_matrix_file, "d_matrix")
        catch e
            @warn "Failed to load cached distance matrix from $d_matrix_file: $e. Recomputing..."
            # Explicit recomputation path if load fails
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

"""
    precompute_dimensions(datadir::String, mutations::Vector{String}, round::Int=1)

Precompute dimensions for a list of mutations.
"""
function precompute_dimensions(datadir::String, mutations::Vector{String}, round_val::Int=1) # Renamed round
    dimensions = Dict{String, Int}()
    for mutation_id in mutations # Changed `mutation` to `mutation_id` for clarity within loop
        try
            dim = get_structure_dimension(datadir, mutation_id, round_val)
            dimensions[mutation_id] = dim
            println("Calculated dimension for $mutation_id: $dim")
        catch e
            @warn "Failed to compute dimension for $mutation_id: $e"
        end
    end
    return dimensions
end

"""
    process_multiple_mutations(datadir::String, mutations::Vector{String}, round::Int; parallel::Bool=false)

Process multiple mutations and calculate distance maps.
"""
function process_multiple_mutations(datadir::String, mutations::Vector{String}, 
                                   round_val::Int; # Renamed round
                                   parallel::Bool=false)
    results = Dict{String, Matrix{Float64}}()
    
    if parallel
        # Ensure Distributed is loaded and workers are available.
        # Functions called by remotecall_fetch must be available on workers.
        # (e.g., via @everywhere include("this_file.jl") or @everywhere using ModuleName)
        num_workers_available = nworkers()
        if num_workers_available == 0
            @warn "Parallel processing requested, but no worker processes available. Falling back to sequential."
            # Fallback handled by subsequent `else` for sequential processing
        else
            # Using original logic for distributing tasks
            tasks = [@spawn remotecall_fetch(get_dist_map, workers()[mod1(i, num_workers_available)], datadir, mutation, round_val) for (i, mutation) in enumerate(mutations)]
            fetched_results = fetch.(tasks) # Collect results
            
            for (i, mutation_id) in enumerate(mutations)
                results[mutation_id] = fetched_results[i]
                # Optional: Add print statement for progress if desired
                # println("Processed $mutation_id in parallel")
            end
            return results # Return after parallel processing
        end
    end
    
    # Sequential processing (also fallback if parallel setup fails)
    for mutation_id in mutations # Changed `mutation` to `mutation_id` for clarity
        results[mutation_id] = get_dist_map(datadir, mutation_id, round_val)
        println("Processed $mutation_id") # Original print statement
    end
    
    return results
end

"""
    get_low_plddt_residues(mutation::String, round::Int, datadir::String; threshold::Float64=90.0)

Identify residues with low pLDDT scores.
"""
function get_low_plddt_residues(mutation::String, round_val::Int, datadir::String; threshold::Float64=90.0) # Renamed round
    subfolder_pattern = "$(datadir)/$(mutation)_$(round_val)/seed-*_sample-0"
    subfolders = glob(subfolder_pattern)
    
    if isempty(subfolders)
        error("No matching folder found: $subfolder_pattern")
    end
    
    subfolder = subfolders[1]
    model_pdb = joinpath(subfolder, "model.pdb")
    # temp_pdb = "$(mutation)_temp_pl.pdb" # This variable was unused in the original.
    
    _ensure_pdb_from_cif_if_needed(subfolder, model_pdb) # Handles PDB creation

    structure = read(model_pdb, PDBFormat)
    
    confidence_file = joinpath(subfolder, "confidences.json")
    if !isfile(confidence_file)
        error("Confidence file not found: $confidence_file")
    end
    
    confidence_data = JSON.parsefile(confidence_file)
    # Ensure atom_plddts is Vector{Float64} for `mean` function compatibility
    atom_plddts = convert(Vector{Float64}, confidence_data["atom_plddts"])
    
    residue_plddts = Dict{Int, Float64}()
    atom_idx_tracker = 1 # Use a different name to avoid confusion with atom_index if it's a global/module var
    
    for model_obj in structure # model is a type in BioStructures
        for chain_obj in model_obj # chain is a type in BioStructures
            for residue_obj in chain_obj # residue is a type in BioStructures
                # Original parsing logic, will error if resid is not parsable to Int
                res_id = parse(Int, BioStructures.resid(residue_obj))
                
                # Using collectatoms for consistency, original had collect(atoms(residue))
                current_residue_atoms = collectatoms(residue_obj)
                atom_count = length(current_residue_atoms)
                
                if atom_count > 0 && !haskey(residue_plddts, res_id)
                    # Ensure atom_idx_tracker and atom_count do not exceed atom_plddts bounds
                    if atom_idx_tracker + atom_count - 1 <= length(atom_plddts)
                        plddt_scores_for_residue = view(atom_plddts, atom_idx_tracker:(atom_idx_tracker + atom_count - 1))
                        residue_plddts[res_id] = mean(plddt_scores_for_residue)
                    else
                        @warn "Atom index out of bounds for residue $res_id (mutation $mutation, round $round_val). Atom pLDDTs length: $(length(atom_plddts)), trying to access up to $(atom_idx_tracker + atom_count - 1). Skipping pLDDT for this residue."
                    end
                elseif atom_count == 0 && !haskey(residue_plddts, res_id) # Handle residues with no atoms if they occur
                     @warn "Residue $res_id (mutation $mutation, round $round_val) has no atoms. Skipping pLDDT calculation."
                end
                
                # Crucially, advance atom_idx_tracker by atom_count for every residue processed,
                # regardless of whether its pLDDT was stored, to maintain sync with atom_plddts list.
                # This fixes a potential bug in the original atom_index advancement.
                atom_idx_tracker += atom_count
            end
        end
    end
    
    low_plddt_res_ids = [res_id for (res_id, plddt) in residue_plddts if plddt < threshold]
    
    return sort(low_plddt_res_ids)
end

"""
    find_residues_within_distance(residue_id::Int, d_matrix::Matrix{Float64}; distance::Float64=13.0)

Find all residues within a specific distance of a target residue.
Original logic: `x -> x < distance`. This includes the residue itself (dist 0)
and potentially others if distance is very small.
"""
function find_residues_within_distance(residue_id::Int, d_matrix::Matrix{Float64}; distance::Float64=13.0)
    if !(1 <= residue_id <= size(d_matrix, 1))
        @warn "Residue ID $residue_id is out of bounds for distance matrix row access (max rows: $(size(d_matrix,1))). Returning empty list."
        return Int[]
    end
    # Original logic: findall(x -> x < distance, d_matrix[residue_id, :])
    # This includes self-distance (0.0) if distance > 0.
    # Subsequent code (calculate_strain, collect_strains) handles `if residue == site; continue; end`.
    nearby_indices = findall(x -> x < distance, d_matrix[residue_id, :])
    
    return sort(nearby_indices) # Returns sorted Vector{Int}
end

"""
    read_pae(datadir::String, mutation::String, round::Int)

Read the Predicted Aligned Error (PAE) matrix.
"""
function read_pae(datadir::String, mutation::String, round_val::Int) # Renamed round
    subfolder_pattern = "$(datadir)/$(mutation)_$(round_val)/seed-*_sample-0"
    subfolders = glob(subfolder_pattern)
    
    if isempty(subfolders)
        error("No matching folder found for $mutation round $round_val")
    end
    
    subfolder = subfolders[1]
    confidence_file = joinpath(subfolder, "confidences.json")
    
    if !isfile(confidence_file)
        error("Confidence file not found for $mutation round $round_val: $confidence_file")
    end
    
    try
        raw_data = JSON.parsefile(confidence_file)
        pae_any = raw_data["pae"] # Expected: Vector of Vectors (rows)
        
        # Original conversion logic using map and hcat
        pae_matrix_rows = map(row -> map(Float64, row), pae_any) # Convert inner elements to Float64
        pae_matrix = Matrix(hcat(pae_matrix_rows...)') # Convert Vector{Vector} to Matrix
        
        return pae_matrix
        
    catch e
        # Original error message structure
        error("Error processing PAE data for $mutation round $round_val: $e")
    end
end

"""
    calculate_strain(alpha::Float64, pae::Matrix{Float64}, d_matrix::Matrix{Float64}, site::Int64)

Standalone function to calculate strain for a single site using given PAE and distance matrix.
This was present in the original file and is kept as is.
The denominator uses length(nearby_residues) which includes 'site' itself if `distance > 0`.
"""
function calculate_strain(alpha::Float64, pae::Matrix{Float64}, d_matrix::Matrix{Float64}, site::Int64)
    # `find_residues_within_distance` uses default distance of 13.0 Angstroms
    # If a different distance is intended here, it needs to be passed.
    # Assuming default is fine as per original structure.
    nearby_residues = find_residues_within_distance(site, d_matrix) 
    strain = 0.0
    
    # Check if pae and d_matrix are large enough for site access
    # Basic check, assumes square matrices for simplicity or symmetric access patterns.
    max_dim_pae = minimum(size(pae))
    max_dim_dists = minimum(size(d_matrix))

    valid_neighbors_for_sum = 0 # To count actual terms added to strain sum
    for residue_neighbor in nearby_residues
        if residue_neighbor == site
            continue
        end
        
        # Bounds check before accessing matrices
        if site > size(pae,1) || residue_neighbor > size(pae,2) || 
           site > size(d_matrix,1) || residue_neighbor > size(d_matrix,2)
            @warn "Site $site or neighbor $residue_neighbor out of bounds. PAE dims: $(size(pae)), Dists dims: $(size(d_matrix)). Skipping pair."
            continue
        end

        dist_val = d_matrix[site, residue_neighbor]
        if dist_val <= 0
            @warn "Non-positive distance ($dist_val) for site $site, residue $residue_neighbor. Skipping."
            continue
        end

        strain += pae[site, residue_neighbor] / (dist_val^alpha)
        valid_neighbors_for_sum += 1
    end

    # Denominator logic based on original `calculate_strain`'s use of `length(nearby_residues)`
    # which counts all residues within sphere (including self if dist_val < distance_cutoff).
    # However, the summation loop skips self.
    # If nearby_residues is empty or only contains `site`, length might be 0 or 1.
    # The original `collect_strains` used a `valid_count` for its average.
    # This standalone `calculate_strain` is a bit ambiguous.
    # Using `valid_neighbors_for_sum` for the average is more robust if that's the intent.
    # If `length(nearby_residues)` must be used:
    # num_denominator = length(nearby_residues)
    # If `site` is in `nearby_residues`, and we want to average over *other* neighbors,
    # then num_denominator should be `length(nearby_residues) - 1` (if site was present).
    # The original `collect_strains` averaging method (sum / N_interacting_partners) is clearer.
    # For this standalone function, sticking to `length(nearby_residues)` as per original structure.
    # Add check for division by zero.
    if length(nearby_residues) > 0
        strain = strain / length(nearby_residues) 
    else
        strain = 0.0 # Or handle as an error, or return NaN
        @warn "No residues found by find_residues_within_distance for site $site to average strain. Strain set to 0."
    end

    return strain
end


"""
    collect_strains(datadir::String, protein::String; ...)

Collect and average strain values for multiple residues across multiple rounds.
"""
function collect_strains(datadir::String, protein::String; 
                         alpha::Float64=2.0, 
                         residue_range=nothing, # Original name
                         num_rounds::Int=20,
                         verbose::Bool=true,
                         cache::Bool=true, # Original name
                         plddt_threshold::Float64=90.0)
    
    effective_residue_range = residue_range # Use a new variable for modifications
    if effective_residue_range === nothing
        try
            dists_r1 = get_dist_map(datadir, protein, 1)
            # Ensure it's a collection (e.g., Vector) for `filter` compatibility later
            effective_residue_range = collect(1:size(dists_r1, 1)) 
            if verbose
                println("Residue range not provided, using full range: $(first(effective_residue_range)) to $(last(effective_residue_range))")
            end
        catch e
            @error "Failed to determine structure dimension for $protein (round 1) to set residue range: $e. Cannot proceed."
            return Dict{Int, Float64}()
        end
    else
        # Ensure it's a collection if a range like 10:20 was passed
        effective_residue_range = collect(effective_residue_range)
    end
    
    if isempty(effective_residue_range)
        @warn "Initial residue range for $protein is empty. No strains will be calculated."
        return Dict{Int, Float64}()
    end

    try
        low_plddt_residues = get_low_plddt_residues(protein, 1, datadir, threshold=plddt_threshold)
        if !isempty(low_plddt_residues)
            original_count = length(effective_residue_range)
            # Filter out low pLDDT residues
            effective_residue_range = filter(r -> !(r in low_plddt_residues), effective_residue_range)
            
            if verbose
                excluded_count = original_count - length(effective_residue_range)
                println("Excluded $excluded_count residues with pLDDT < $plddt_threshold (based on round 1 of $protein)")
            end
        end
    catch e
        @warn "Failed to filter low pLDDT residues for $protein (round 1): $e. Proceeding with current residue list."
    end

    if isempty(effective_residue_range)
        @warn "Residue list for $protein became empty after pLDDT filtering. No strains will be calculated."
        return Dict{Int, Float64}()
    end
    
    cache_dir = joinpath(datadir, "strain_cache")
    if cache && !isdir(cache_dir)
        try mkpath(cache_dir) catch e
            @warn "Failed to create cache directory $cache_dir: $e. Disabling cache for this run."
            cache = false # Update local cache flag
        end
    end
    
    site_value_sums = Dict(k => 0.0 for k in effective_residue_range)
    
    # This cache stores neighbors based on the structure of the round where the site is first encountered.
    nearby_residues_cache = Dict{Int, Vector{Int}}() 
    
    processed_rounds_count = 0 # To average correctly if some rounds fail entirely

    # --- Inner function to process one round, mimicking original structure ---
    function process_round_logic(r_idx::Int) # Renamed 'r' to avoid conflict
        local_results_for_round = Dict(k => 0.0 for k in effective_residue_range)
        
        plddt_suffix = "_plddt$(plddt_threshold)" # Original cache suffix
        # Alpha is a float, format it consistently for filename
        alpha_str = replace(string(alpha), "." => "p") 
        cache_file = joinpath(cache_dir, "$(protein)_round$(r_idx)_alpha$(alpha_str)$(plddt_suffix).jld2")
        
        if cache && isfile(cache_file)
            if verbose; println("Loading cached data for $protein round $r_idx from $cache_file"); end
            try
                # Load and ensure all keys in effective_residue_range are populated
                loaded_data = load(cache_file, "strain_data")
                for site_key in effective_residue_range
                    local_results_for_round[site_key] = get(loaded_data, site_key, 0.0) # Default if missing
                end
                return local_results_for_round, true # Indicate success
            catch e
                @warn "Failed to load cached data for $protein round $r_idx from $cache_file: $e. Recomputing..."
            end
        end
        
        try
            dists_current_round = get_dist_map(datadir, protein, r_idx)
            paes_current_round = read_pae(datadir, protein, r_idx)
            
            for site_id in effective_residue_range
                # Populate nearby_residues_cache if site_id is new.
                # Uses dists_current_round, so neighbors are from this round's structure if new.
                if !haskey(nearby_residues_cache, site_id)
                    # Default distance for find_residues_within_distance is 13.0.
                    # If a different cutoff is desired for defining neighbors, it should be a parameter.
                    nearby_residues_cache[site_id] = find_residues_within_distance(site_id, dists_current_round)
                end
                
                current_site_neighbors = nearby_residues_cache[site_id]
                
                strain_sum_for_site = 0.0
                valid_interacting_partners = 0
                
                for neighbor_id in current_site_neighbors
                    if neighbor_id == site_id
                        continue 
                    end
                    
                    # Bounds checking for matrix access
                    if !(1 <= site_id <= size(paes_current_round, 1) && 1 <= neighbor_id <= size(paes_current_round, 2) &&
                         1 <= site_id <= size(dists_current_round, 1) && 1 <= neighbor_id <= size(dists_current_round, 2))
                        @warn "Site $site_id or neighbor $neighbor_id out of bounds for PAE/distance matrices in $protein round $r_idx. Skipping pair."
                        continue
                    end
                    
                    distance_val = dists_current_round[site_id, neighbor_id]
                    if distance_val <= 0.0 # Should not happen for distinct CAs
                        @warn "Non-positive distance ($distance_val) between site $site_id and $neighbor_id in $protein round $r_idx. Skipping pair."
                        continue
                    end

                    strain_sum_for_site += paes_current_round[site_id, neighbor_id] / (distance_val^alpha)
                    valid_interacting_partners += 1
                end
                
                if valid_interacting_partners > 0
                    local_results_for_round[site_id] = strain_sum_for_site / valid_interacting_partners
                else
                    local_results_for_round[site_id] = 0.0 # Keep as 0.0 if no valid interactions
                    if !isempty(current_site_neighbors) # Warn only if neighbors existed but were invalid
                        @warn "No valid interacting partners found for site $site_id in $protein round $r_idx, though neighbors were identified. Strain set to 0."
                    end
                end
            end
            
            if cache
                try
                    jldsave(cache_file; strain_data=local_results_for_round)
                catch e
                    @warn "Failed to cache results for $protein round $r_idx to $cache_file: $e"
                end
            end
            
            if verbose; println("Completed processing for $protein round $r_idx."); end
            return local_results_for_round, true # Indicate success

        catch e
            @error "Failed to process $protein round $r_idx: $e. Strains for this round will be 0.0 for all sites."
            # Return dict of zeros for this round on error
            return Dict(k => 0.0 for k in effective_residue_range), false # Indicate failure
        end
    end # --- End of process_round_logic ---
    
    # Process rounds sequentially and aggregate
    for r_val in 1:num_rounds # Renamed r to r_val
        if verbose; println("Processing $protein, round $r_val of $num_rounds..."); end
        
        round_strains, success = process_round_logic(r_val)
        
        if success
            processed_rounds_count += 1
            for (site, strain_val) in round_strains
                if haskey(site_value_sums, site) # Ensure site is in our target list
                    site_value_sums[site] += strain_val
                end
            end
        end
    end
    
    final_averaged_values = Dict{eltype(effective_residue_range), Float64}()
    if processed_rounds_count > 0
        for site_id_key in effective_residue_range # Iterate over the keys we expect
            final_averaged_values[site_id_key] = site_value_sums[site_id_key] / processed_rounds_count
        end
    else
        @warn "No rounds were successfully processed for $protein. Averaged strains will be 0.0 or based on initial sums."
        # Populate with zeros if no rounds processed, or sums remain 0.0
        for site_id_key in effective_residue_range
            final_averaged_values[site_id_key] = 0.0
        end
    end
        
    if verbose
        println("\nProcessing completed for $protein.")
        println("Averaged over $processed_rounds_count successfully processed rounds (out of $num_rounds attempted).")
        println("\nFinal Averaged Strain Values per Site:")
        # Sort by site ID for consistent output
        sorted_sites = sort(collect(keys(final_averaged_values)))
        for site_print in sorted_sites
            println("Site $site_print: $(final_averaged_values[site_print])")
        end
    end
    
    return final_averaged_values
end

"""
    calculate_ME(Strain::Dict{Int, Float64}, Strain_wt::Dict{Int, Float64})

Calculate Mutation Entropy (ME).
"""
function calculate_ME(Strain::Dict{Int, Float64}, Strain_wt::Dict{Int, Float64})
    common_sites = intersect(keys(Strain), keys(Strain_wt))
    
    if isempty(common_sites)
        # Original behavior: error.
        error("No common residue sites found between the two strain dictionaries")
    end
    
    ME_vals = Dict{Int, Float64}() # Renamed ME to ME_vals for clarity
    for site_id in common_sites # Renamed site to site_id
        ME_vals[site_id] = abs(Strain[site_id] - Strain_wt[site_id])
    end
    
    # Original logging structure
    num_strain_sites = length(keys(Strain))
    num_strain_wt_sites = length(keys(Strain_wt))
    num_common_sites = length(common_sites)
    
    missing_from_ME = (num_strain_sites - num_common_sites) + (num_strain_wt_sites - num_common_sites)
    if missing_from_ME > 0 # If there were sites not in common
        @info "$(missing_from_ME) residue sites were excluded from ME calculation due to missing in either strain dataset"
        @info "Calculated ME for $(num_common_sites) common residue sites"
    else
         @info "Calculated ME for $(num_common_sites) common residue sites (all sites matched)." # Added for completeness
    end
    
    return ME_vals
end