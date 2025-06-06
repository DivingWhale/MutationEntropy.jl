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
    nearby_residues = MutationEntropy.find_residues_within_distance(mutation_site, dists, distance=13.0)
    
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
    
    # Highlight regions with residues within 13Å of the mutation site
    for residue in residues
        if residue in nearby_residues
            # Draw a vertical band to highlight these residues
            vspan!(ax, residue-0.5, residue+0.5, color=("#FFBD69", 0.3))
        end
    end
    
    # Line plot of ME values
    lines!(ax, residues, me_values, 
           color="#381460", 
           linewidth=2)
    
    # Add scatter points for each residue
    scatter!(ax, residues, me_values, 
             color="#B21F66", 
             markersize=5)
    
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
    lines!(ax, [NaN], [NaN], color=:transparent, label="Residues within 13Å")
    poly!(ax, [NaN, NaN, NaN], [NaN, NaN, NaN], color=("#FFBD69", 0.3), strokewidth=0)
    
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

# Main function to calculate and plot ME for a specific variant
function main()
    # Parameters
    data_path = "data"
    mutant = "m147g"  # Example mutant
    alpha = 3.0       # Alpha parameter for strain calculation
    residue_range = 83:231  # Range of residues to analyze
    num_rounds = 20   # Number of rounds to average over
    
    println("Starting ME calculation for $mutant with residue range $residue_range")
    
    # Calculate wild-type strain values
    println("Calculating wild-type strain values...")
    wt_avg_S_values = MutationEntropy.collect_strains(
        data_path,                # Base directory for data
        "thermonulease",          # Protein name
        alpha=alpha,              # Alpha parameter for strain calculation
        residue_range=residue_range, # Range of residues to analyze
        num_rounds=num_rounds,    # Number of rounds to average
        verbose=false,             # Print progress information
        cache=true,               # Use cache if available
        plddt_threshold=95.0      # pLDDT threshold below which residues are excluded
    )
    
    # Calculate mutant strain values
    println("\nCalculating mutant strain values for $mutant...")
    mutant_avg_S_values = MutationEntropy.collect_strains(
        data_path,
        mutant,
        alpha=alpha,
        residue_range=residue_range,
        num_rounds=num_rounds,
        verbose=false,
        cache=true,
        plddt_threshold=95.0
    )
    
    # Calculate Mutation Entropy
    println("\nCalculating ME...")
    mutant_ME = MutationEntropy.calculate_ME(mutant_avg_S_values, wt_avg_S_values)
    
    # Get statistics
    me_values = collect(values(mutant_ME))
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
    
    # Create ME vs Residue plot
    println("\nCreating ME vs Residue plot...")
    fig = plot_ME_vs_residue(mutant_ME, mutant, residue_range=residue_range)
    
    # Also create standard ME vs Distance plot for comparison
    println("\nCreating standard ME vs Distance plot...")
    dist_fig = MutationEntropy.plot_MEvsDist(mutant_ME, data_path, mutant)
    save_path = joinpath("figs", "ME_vs_Dist_$(mutant).png")
    save(save_path, dist_fig)
    println("Distance plot saved to: $save_path")
    
    println("\nAnalysis complete!")
end

# Run the main function
main()
