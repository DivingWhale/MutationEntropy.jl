"""
    ca_distMap(struc)

Calculate distance map between C-alpha atoms in the given protein structure.
"""
function ca_distMap(struc) 
    dists = DistanceMap(collectatoms(struc, calphaselector))
    return dists
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

    if length(nearby_residues) > 1
        strain = strain / (length(nearby_residues) - 1) # Exclude self from average
    else
        strain = 0.0 # Or handle as an error, or return NaN
        @warn "No valid neighbors found for site $site to average strain. Strain set to 0."
    end

    return strain
end


"""
    collect_strains(datadir::String, protein::String; ...)

Collect and average strain values for multiple residues across multiple rounds.
"""
function collect_strains(datadir::String, protein::String; 
                         alpha::Float64=2.0, 
                         residue_range=nothing, 
                         num_rounds::Int=20,
                         verbose::Bool=true,
                         cache::Bool=true, 
                         plddt_threshold::Float64=90.0)
    
    effective_residue_range = residue_range
    if effective_residue_range === nothing
        try
            dists_r1 = get_dist_map(datadir, protein, 1)
            effective_residue_range = collect(axes(dists_r1, 1))
            verbose && println("Residue range not provided, using full range: $(first(effective_residue_range)) to $(last(effective_residue_range))")
        catch e
            @error "Failed to determine structure dimension for $protein (round 1) to set residue range: $e. Cannot proceed."
            return Dict{Int, Float64}()
        end
    else
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
            cache = false
        end
    end
    
    site_value_sums = Dict(k => 0.0 for k in effective_residue_range)
    processed_rounds_count = 0

    for r_idx in 1:num_rounds
        verbose && println("Processing $protein, round $r_idx of $num_rounds...")
        
        plddt_suffix = "_plddt$(plddt_threshold)"
        alpha_str = replace(string(alpha), "." => "p") 
        cache_file = joinpath(cache_dir, "$(protein)_round$(r_idx)_alpha$(alpha_str)$(plddt_suffix).jld2")
        
        round_strains = Dict{Int, Float64}()
        if cache && isfile(cache_file)
            verbose && println("Loading cached data for $protein round $r_idx from $cache_file")
            try
                loaded_data = load(cache_file, "strain_data")
                for site_key in effective_residue_range
                    round_strains[site_key] = get(loaded_data, site_key, 0.0)
                end
            catch e
                @warn "Failed to load cached data for $protein round $r_idx from $cache_file: $e. Recomputing..."
            end
        end

        if isempty(round_strains)
            try
                dists_current_round = get_dist_map(datadir, protein, r_idx)
                paes_current_round = read_round_pae(datadir, protein, r_idx)
                
                for site_id in effective_residue_range
                    round_strains[site_id] = calculate_strain(alpha, paes_current_round, dists_current_round, site_id)
                end
                
                if cache
                    try
                        jldsave(cache_file; strain_data=round_strains)
                    catch e
                        @warn "Failed to cache results for $protein round $r_idx to $cache_file: $e"
                    end
                end
                verbose && println("Completed processing for $protein round $r_idx.")
            catch e
                @error "Failed to process $protein round $r_idx: $e. Strains for this round will be 0.0 for all sites."
                for site_id in effective_residue_range
                    round_strains[site_id] = 0.0
                end
            end
        end
        
        processed_rounds_count += 1
        for (site, strain_val) in round_strains
            if haskey(site_value_sums, site)
                site_value_sums[site] += strain_val
            end
        end
    end
    
    final_averaged_values = Dict{eltype(effective_residue_range), Float64}()
    if processed_rounds_count > 0
        for site_id_key in effective_residue_range
            final_averaged_values[site_id_key] = site_value_sums[site_id_key] / processed_rounds_count
        end
    else
        @warn "No rounds were successfully processed for $protein. Averaged strains will be 0.0 or based on initial sums."
        for site_id_key in effective_residue_range
            final_averaged_values[site_id_key] = 0.0
        end
    end
        
    if verbose
        println("\nProcessing completed for $protein.")
        println("Averaged over $processed_rounds_count successfully processed rounds (out of $num_rounds attempted).")
        println("\nFinal Averaged Strain Values per Site:")
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