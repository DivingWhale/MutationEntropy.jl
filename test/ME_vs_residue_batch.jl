using CairoMakie
using MutationEntropy
using Statistics

"""
    plot_ME_vs_residue(ME::Dict{Int, Float64}, mutation::String; residue_range=nothing)

Plot Mutation Entropy (ME) against residue numbers as a line graph.
Adds a vertical line at the mutation position.

# Arguments
- `ME::Dict{Int, Float64}`: Dictionary mapping residue IDs to their mutation entropy values
- `mutation::String`: Mutation identifier (e.g., "a140e")
- `residue_range`: Optional range of residues to include in the plot (default: all residues in ME)

# Returns
- The figure object
"""
function plot_ME_vs_residue(ME::Dict{Int, Float64}, mutation::String; residue_range=nothing)
    # Extract mutation site residue ID (e.g., get 140 from "a140e")
    mutation_site = parse(Int, match(r"[a-zA-Z](\d+)[a-zA-Z]", mutation).captures[1])
    
    # Extract residue IDs and mutation entropy values
    if residue_range === nothing
        residues = sort(collect(keys(ME)))
    else
        residues = sort(filter(r -> r in keys(ME) && r in residue_range, collect(residue_range)))
    end
    
    me_values = [get(ME, r, NaN) for r in residues]
    
    # Get distance matrix to find residues within 13Å of mutation site
    data_path = "data"
    dists = MutationEntropy.get_dist_map(data_path, mutation, 1)
    max_distance = 13.0
    nearby_residues = MutationEntropy.find_residues_within_distance(mutation_site, dists, distance=max_distance)
    
    # Create a dictionary mapping residue to its distance from mutation site
    residue_distances = Dict{Int, Float64}()
    for residue in nearby_residues
        if residue != mutation_site  # Skip the mutation site itself
            residue_distances[residue] = dists[mutation_site, residue]
        end
    end
    
    # Get list of residues with low pLDDT scores that were excluded
    plddt_threshold = 90.0
    low_plddt_residues = MutationEntropy.get_low_plddt_residues(mutation, 1, data_path, threshold=plddt_threshold)
    # Filter to only include those in our range of interest
    if residue_range !== nothing
        low_plddt_residues = filter(r -> r in residue_range, low_plddt_residues)
    end
    
    # Create a figure with CairoMakie
    fig = Figure(size=(1200, 800))
    ax = Axis(fig[1, 1], 
              title="Mutation Effect vs Residue Number for $(uppercase(mutation))",
              xlabel="Residue Number",
              ylabel="Mutation Effect (ME)")
    
    # Set custom x-axis ticks with intervals of 20
    if !isempty(residues)
        min_res = minimum(residues)
        max_res = maximum(residues)
        start_tick = ceil(min_res/20) * 20  # Round up to nearest multiple of 20
        xticks = start_tick:20:max_res
        ax.xticks = xticks
    end
    
    # Create full range including low pLDDT residues (for plotting excluded points)
    all_residues_in_range = residue_range === nothing ? (minimum(residues):maximum(residues)) : residue_range
    
    # Define distance range for opacity scaling (ignore the mutation site itself)
    min_dist = minimum(values(residue_distances); init=max_distance)
    
    # Define opacity range (0.2 to 0.8)
    min_opacity = 0.2
    max_opacity = 0.8
    
    # Highlight regions with residues within 13Å of the mutation site
    for residue in residues
        if residue in keys(residue_distances)
            # Get the distance and calculate opacity based on distance
            dist = residue_distances[residue]
            # Normalize distance to opacity range (closer = more opaque)
            opacity = max_opacity - (dist - min_dist) / (max_distance - min_dist) * (max_opacity - min_opacity)
            opacity = clamp(opacity, min_opacity, max_opacity)
            
            # Draw a vertical band to highlight these residues with distance-dependent opacity
            vspan!(ax, residue-0.5, residue+0.5, color=("#F97A00", opacity))
        end
    end
    vspan!(ax, mutation_site-0.5, mutation_site+0.5, color=("#F97A00", 0.8))
    
    # Line plot of ME values
    lines!(ax, residues, me_values, 
           color="#381460", 
           linewidth=2)
    
    # Add scatter points for each residue
    scatter!(ax, residues, me_values, 
             color="#381460", 
             markersize=8)
    
    # Add a vertical line at the mutation site
    vlines!(ax, mutation_site, 
            color=:gray, 
            linestyle=:dash, 
            linewidth=2,
            label="Mutation Site ($(mutation_site))")
    
    # Plot excluded points (low pLDDT residues) if any
    if !isempty(low_plddt_residues)
        # Create dummy points at the bottom of the plot for excluded residues
        y_min = minimum(me_values) * 0.9  # Place slightly below the lowest ME value
        excluded_x = filter(r -> r in low_plddt_residues, collect(all_residues_in_range))
        excluded_y = fill(y_min, length(excluded_x))
        
        # Add points for excluded residues
        scatter!(ax, excluded_x, excluded_y, 
                 color="#B21F66", 
                 marker=:xcross,
                 markersize=8,
                 label="Excluded (pLDDT < 90)")
    end
    
    # Add a text label for the mutation site
    text!(ax, mutation_site, maximum(me_values) * 1.05, 
          text="$(uppercase(mutation))", 
          align=(:center, :bottom), 
          fontsize=12)
    
    # Add legend with labels for highlighted region
    lines!(ax, [NaN], [NaN], color=:transparent, label="Residues within 13Å (darker = closer)")
    # Create a gradient rectangle for the legend to represent distance
    poly!(ax, [NaN, NaN, NaN], [NaN, NaN, NaN], color=("#FFBD69", 0.5), strokewidth=0)
    
    # Add legend - use a valid position symbol that works with CairoMakie
    axislegend(ax, position=:rt)  # :rt means right-top
    
    # Save the plot
    save_dir = "figs"
    if !isdir(save_dir)
        mkpath(save_dir)
    end
    save_path = joinpath(save_dir, "ME_vs_residue_$(mutation).png")
    save(save_path, fig)
    println("Plot saved to: $save_path")
    
    return fig
end

"""
    process_mutant(mutant::String, wt_avg_S_values::Dict{Int, Float64}; ...)
    
Process a single mutant and generate plots for Mutation Entropy.

# Arguments
- `mutant::String`: Mutant identifier (e.g., "a140e")
- `wt_avg_S_values::Dict{Int, Float64}`: Wild-type strain values
- `data_path::String="data"`: Base directory containing protein data
- `alpha::Float64=3.0`: Alpha parameter for strain calculation
- `residue_range=83:231`: Range of residues to analyze
- `num_rounds::Int=20`: Number of rounds to average over
- `verbose::Bool=true`: Whether to print progress information
- `plddt_threshold::Float64=95.0`: pLDDT threshold below which residues are excluded

# Returns
- Dictionary mapping residue IDs to their mutation entropy values
"""
function process_mutant(mutant::String, wt_avg_S_values::Dict{Int, Float64}; 
                       data_path::String="data", 
                       alpha::Float64=3.0, 
                       residue_range=83:231, 
                       num_rounds::Int=20, 
                       verbose::Bool=true,
                       plddt_threshold::Float64=95.0)
    
    if verbose
        println("\n" * "="^50)
        println("Processing mutant: $mutant")
        println("="^50)
        println("Calculating mutant strain values...")
    end
    
    # Calculate mutant strain values
    mutant_avg_S_values = MutationEntropy.collect_strains(
        data_path,
        mutant,
        alpha=alpha,
        residue_range=residue_range,
        num_rounds=num_rounds,
        verbose=verbose,
        cache=true,
        plddt_threshold=plddt_threshold
    )
    
    # Calculate Mutation Entropy
    if verbose
        println("\nCalculating Mutation Entropy (ME)...")
    end
    mutant_ME = MutationEntropy.calculate_ME(mutant_avg_S_values, wt_avg_S_values)
    
    # Get statistics
    me_values = collect(values(mutant_ME))
    if verbose
        println("ME Statistics:")
        println("  - Number of residues with ME: $(length(mutant_ME))")
        println("  - Mean ME: $(mean(me_values))")
        
        # Find residue with max ME
        max_me_val = maximum(me_values)
        max_me_residue = [k for (k,v) in mutant_ME if v == max_me_val][1]
        println("  - Max ME: $(max_me_val) at residue $(max_me_residue)")
        
        # Find residue with min ME
        min_me_val = minimum(me_values)
        min_me_residue = [k for (k,v) in mutant_ME if v == min_me_val][1]
        println("  - Min ME: $(min_me_val) at residue $(min_me_residue)")
    end
    
    # Create ME vs Residue plot
    if verbose
        println("\nCreating ME vs Residue plot...")
    end
    fig = plot_ME_vs_residue(mutant_ME, mutant, residue_range=residue_range)
    
    # Also create standard ME vs Distance plot for comparison
    if verbose
        println("\nCreating standard ME vs Distance plot...")
    end
    dist_fig = MutationEntropy.plot_MEvsDist(mutant_ME, data_path, mutant)
    save_path = joinpath("figs", "ME_vs_Dist_$(mutant).png")
    save(save_path, dist_fig)
    if verbose
        println("Distance plot saved to: $save_path")
    end
    
    return mutant_ME
end

# Main function to calculate and plot ME for multiple variants
function main()
    # Parameters
    data_path = "data"
    mutants = ["a140e", "l120e", "e155f", "e183f", "f143g", "m147g", "i154g"] # List of mutants to process
    alpha = 3.0       # Alpha parameter for strain calculation
    residue_range = 83:231  # Range of residues to analyze
    num_rounds = 20   # Number of rounds to average over
    plddt_threshold = 90.0  # pLDDT threshold below which residues are excluded
    
    println("Starting ME calculation for $(length(mutants)) mutants with residue range $residue_range")
    println("Alpha = $alpha, pLDDT threshold = $plddt_threshold")
    
    # Create figs directory if it doesn't exist
    if !isdir("figs")
        mkpath("figs")
    end
    
    # Calculate wild-type strain values (only need to do this once)
    println("\nCalculating wild-type strain values...")
    wt_avg_S_values = MutationEntropy.collect_strains(
        data_path,                # Base directory for data
        "thermonulease",          # Protein name
        alpha=alpha,              # Alpha parameter for strain calculation
        residue_range=residue_range, # Range of residues to analyze
        num_rounds=num_rounds,    # Number of rounds to average
        verbose=false,             # Print progress information
        cache=true,               # Use cache if available
        plddt_threshold=plddt_threshold
    )
    
    # Process each mutant
    results = Dict{String, Dict{Int, Float64}}()
    
    for (i, mutant) in enumerate(mutants)
        println("\nProcessing mutant $i/$(length(mutants)): $mutant")
        
        results[mutant] = process_mutant(
            mutant, 
            wt_avg_S_values;
            data_path=data_path,
            alpha=alpha,
            residue_range=residue_range,
            num_rounds=num_rounds,
            verbose=false,
            plddt_threshold=plddt_threshold
        )
    end
    
    println("\nAnalysis complete! Processed $(length(mutants)) mutants:")
    for mutant in mutants
        println("  - $mutant")
    end
    
    return results
end

# Run the main function
main()
