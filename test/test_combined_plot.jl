using CairoMakie
using MutationEntropy
using Printf
using Colors
using Random
using Statistics
using StatsBase

"""
    calculate_combined_ME_data(mutations, data_path, alpha, residue_range, num_rounds, verbose)

Calculate Mutation Entropy (ME) data for multiple protein mutants and save the data to a file.

# Arguments
- `mutations::Vector{String}`: List of mutation identifiers to analyze
- `data_path::String="data"`: Base directory containing protein data
- `alpha::Float64=2.0`: Exponent parameter for strain calculation
- `residue_range=83:231`: Range of residues to analyze
- `num_rounds::Int=20`: Number of simulation rounds to average over
- `verbose::Bool=true`: Whether to print progress information
- `save_data::Bool=true`: Whether to save the calculated data to a file

# Returns
- `Tuple` containing all_distances and all_me_values
"""
function calculate_combined_ME_data(
    mutations::Vector{String}, 
    data_path::String="data", 
    alpha::Float64=2.0; 
    residue_range=83:231, 
    num_rounds::Int=20, 
    verbose::Bool=true,
    save_data::Bool=true
)
    # Calculate wild-type strain values (only need to do this once as reference)
    println("Calculating wild-type strain values...")
    wt_avg_S_values = MutationEntropy.collect_strains(
        data_path,                # Base data directory
        "thermonulease",          # Protein name
        alpha=alpha,              # Alpha parameter for strain calculation
        residue_range=residue_range, # Range of residues to analyze
        num_rounds=num_rounds,    # Number of rounds to average
        verbose=verbose           # Show progress information
    )
    
    # Storage for all data points to be collected
    all_distances = Float64[]
    all_me_values = Float64[]
    mutation_indices = Int[]  # To track which mutation each point belongs to
    
    # Process each mutant
    for (i, mutant) in enumerate(mutations)
        if i % 10 == 0 || i == 1 || i == length(mutations)
            println("Processing mutant $i/$(length(mutations)): $mutant")
        end
        
        # Calculate mutant strain values
        mutant_avg_S_values = MutationEntropy.collect_strains(
            data_path,            # Base data directory
            mutant,               # Mutant name
            alpha=alpha,          # Alpha parameter for strain calculation
            residue_range=residue_range, # Range of residues to analyze
            num_rounds=num_rounds, # Number of rounds to average
            verbose=false         # Less verbose for many mutations
        )
        
        # Calculate Mutation Entropy
        mutant_ME = MutationEntropy.calculate_ME(mutant_avg_S_values, wt_avg_S_values)
        
        # Extract mutation site residue ID (e.g., get 140 from "A140E")
        mutation_site = parse(Int, match(r"[a-zA-Z](\d+)[a-zA-Z]", mutant).captures[1])
        
        # Get low pLDDT residues (pLDDT < 90.0)
        low_plddt_residues = MutationEntropy.get_low_plddt_residues(mutant, 1, data_path)
        
        # Load cached distance matrix for the mutation
        dists = MutationEntropy.get_dist_map(data_path, mutant, 1)
        
        # Extract residue IDs and mutation entropy values (filter out low pLDDT residues)
        residues = collect(keys(mutant_ME))
        residues = filter(r -> !(r in low_plddt_residues), residues)
        me_values = [mutant_ME[r] for r in residues]
        
        # Get distances from each residue to the mutation site
        distances = [dists[mutation_site, res] for res in residues]
        
        # Store the data for later plotting
        append!(all_distances, distances)
        append!(all_me_values, me_values)
        append!(mutation_indices, fill(i, length(distances)))  # Track mutation index
    end
    
    # Save data to file if requested
    if save_data
        mkpath("figs")
        data_file = joinpath("figs", "ME_data_alpha_$(alpha).csv")
        open(data_file, "w") do io
            # Write header
            println(io, "distance,me_value,mutation_index,mutation_name")
            
            # Write data
            for i in 1:length(all_distances)
                mutation_idx = mutation_indices[i]
                mutation_name = mutations[mutation_idx]
                println(io, "$(all_distances[i]),$(all_me_values[i]),$mutation_idx,$(mutation_name)")
            end
        end
        println("Data saved to: $data_file")
    end
    
    # Calculate some statistics about the data
    println("\nStatistics for combined ME data:")
    println("Total number of data points: $(length(all_me_values))")
    println("Total number of mutations: $(length(mutations))")
    println("Mean ME value: $(@sprintf("%.4f", mean(all_me_values)))")
    println("Max ME value: $(@sprintf("%.4f", maximum(all_me_values)))")
    println("Min ME value: $(@sprintf("%.4f", minimum(all_me_values)))")
    println("Correlation between distance and ME: $(@sprintf("%.4f", cor(all_distances, all_me_values)))")
    
    return all_distances, all_me_values
end

"""
    plot_ME_density_heatmap(all_distances, all_me_values, mutations, alpha)
    
Create a density heatmap from the calculated ME data.

# Arguments
- `all_distances::Vector{Float64}`: Vector of distances from mutation sites
- `all_me_values::Vector{Float64}`: Vector of mutation entropy values
- `mutations::Vector{String}`: List of mutation identifiers
- `alpha::Float64=1.0`: Alpha parameter used in ME calculation
- `y_limit::Union{Float64, Nothing}=nothing`: Upper limit for y-axis, automatically determined if nothing

# Returns
- The figure object
"""
function plot_ME_density_heatmap(
    all_distances::Vector{Float64}, 
    all_me_values::Vector{Float64},
    mutations::Vector{String}=String[];
    alpha::Float64=1.0,
    y_limit::Union{Float64, Nothing}=nothing
)
    # Create a figure to hold the data
    fig = Figure(size=(1500, 800))
    ax = Axis(fig[1, 1],
              title="Mutation Effect vs Distance for Multiple Variants (α=$alpha, n=$(length(mutations)))",
              xlabel="Distance from mutation site (Å)",
              ylabel="Mutation Effect",
              aspect=1.5)
              
    # Convert scatter plot to a density heatmap
    println("Creating density heatmap from $(length(all_distances)) data points...")
    
    # Calculate y_limit early so we can adjust bins accordingly
    actual_y_limit = isnothing(y_limit) ? maximum(all_me_values) : y_limit
    
    # Define the binning for the heatmap
    n_bins_x = min(100, div(length(unique(all_distances)), 5) + 1)
    n_bins_y = min(100, div(length(unique(all_me_values)), 5) + 1)
    println("Using $n_bins_x x bins and $n_bins_y y bins for heatmap")
    
    # Calculate the density using a histogram2d approach
    h = fit(Histogram, (all_distances, all_me_values), 
            nbins=(n_bins_x, n_bins_y))
    
    # Get the edges of the bins
    x_edges = h.edges[1]
    y_edges = h.edges[2]
    
    # Get the midpoints of the bins
    x_mids = [(x_edges[i] + x_edges[i+1])/2 for i in 1:length(x_edges)-1]
    y_mids = [(y_edges[i] + y_edges[i+1])/2 for i in 1:length(y_edges)-1]
    
    # Apply log transformation for better visualization (add small constant to avoid log(0))
    density_data = log10.(h.weights .+ 1)
    
    # Plot the heatmap
    heatmap!(ax, x_mids, y_mids, density_data, 
              colormap=:viridis, 
              colorrange=(0, maximum(density_data)))
    
    # Add a colorbar to show the density scale
    cbar = Colorbar(fig[1, 2], colorrange=(0, maximum(density_data)), 
                    colormap=:viridis,
                    label="log10(count + 1)")
    
    # Create figs directory if it doesn't exist
    mkpath("figs")
    
    # Set axis limits for better visualization
    # Note: actual_y_limit is now calculated earlier in the function
    ylims!(ax, 0, actual_y_limit)  # Start y-axis at 0, with customizable upper limit
    
    # Save the combined plot at higher resolution
    save_path = joinpath("figs", "combined_ME_density_alpha_$(alpha).png")
    save(save_path, fig, px_per_unit=2)  # Higher resolution for detailed plots
    println("Combined density plot saved to: $save_path")
    
    return fig
end

"""
    create_combined_ME_plot(mutations, data_path, alpha; kwargs...)

Combined function that both calculates and plots ME data. This is maintained
for backward compatibility and convenience.

# Arguments
See documentation for calculate_combined_ME_data and plot_ME_density_heatmap.
"""
function create_combined_ME_plot(
    mutations::Vector{String}, 
    data_path::String="data", 
    alpha::Float64=2.0; 
    residue_range=83:231, 
    num_rounds::Int=20, 
    verbose::Bool=true,
    y_limit::Union{Float64, Nothing}=nothing,
    save_data::Bool=true
)
    # Calculate the data
    all_distances, all_me_values = calculate_combined_ME_data(
        mutations, data_path, alpha; 
        residue_range=residue_range, 
        num_rounds=num_rounds, 
        verbose=verbose,
        save_data=save_data
    )
    
    # Plot the data
    return plot_ME_density_heatmap(
        all_distances, all_me_values, mutations;
        alpha=alpha,
        y_limit=y_limit
    )
end

# Example usage
function run_combined_ME_analysis()
    # Define the list of mutations to test - this can be replaced with a much larger list
    mutations = ["a140e", "a140g", "a140v", "a142c", "a142f", "a142g", "a142v", "a151g", "a151t", "a151v"]
    
    # Set alpha value
    alpha = 1.0
    
    # Create the combined plot with random colors
    create_combined_ME_plot(mutations, "data", alpha; use_random_colors=true)
end

"""
    find_all_mutations(data_path::String)

Automatically find all mutation folders in the data directory based on folder naming pattern.
"""
function find_all_mutations(data_path::String)
    # Get all directories in the data path
    all_dirs = readdir(data_path, join=true)
    
    # Filter to only include directories
    all_dirs = filter(isdir, all_dirs)
    
    # Extract mutation names based on directory pattern (e.g., a140g_1, a142v_20)
    mutation_dirs = Dict{String, Int}()
    
    for dir in all_dirs
        # Match pattern like "a140g_1" to extract "a140g"
        m = match(r"([a-z]\d+[a-z])_\d+$", basename(dir))
        if m !== nothing
            mutation = m.captures[1]
            # Count occurrences to ensure we have enough rounds
            if haskey(mutation_dirs, mutation)
                mutation_dirs[mutation] += 1
            else
                mutation_dirs[mutation] = 1
            end
        end
    end
    
    # Filter mutations that have at least 10 rounds available
    valid_mutations = [mutation for (mutation, count) in mutation_dirs if count >= 10]
    
    println("Found $(length(valid_mutations)) mutations with at least 10 rounds each.")
    return valid_mutations
end

"""
    run_full_dataset_analysis(data_path::String="data", alpha::Float64=1.0; plot_only::Bool=false)

Run the analysis on all available mutations in the data directory.
This is suitable for large datasets with many variants.

# Arguments
- `data_path::String="data"`: Base directory containing protein data
- `alpha::Float64=1.0`: Alpha parameter for strain calculation
- `plot_only::Bool=false`: Whether to skip calculation and just plot from saved data
- `y_limit::Union{Float64, Nothing}=nothing`: Upper limit for y-axis, automatically determined if nothing
"""
function run_full_dataset_analysis(
    data_path::String="data", 
    alpha::Float64=1.0; 
    plot_only::Bool=false,
    y_limit::Union{Float64, Nothing}=nothing
)
    # Automatically find all mutations in the data directory
    all_mutations = find_all_mutations(data_path)
    
    println("Found $(length(all_mutations)) mutations to analyze.")
    if length(all_mutations) > 0
        if plot_only
            # Try to load saved data
            data_file = joinpath("figs", "ME_data_alpha_$(alpha).csv")
            if isfile(data_file)
                println("Loading data from $data_file...")
                
                # Load the data
                all_distances = Float64[]
                all_me_values = Float64[]
                
                open(data_file, "r") do io
                    # Skip header
                    header = readline(io)
                    
                    # Read data
                    for line in eachline(io)
                        fields = split(line, ",")
                        dist = parse(Float64, fields[1])
                        me_val = parse(Float64, fields[2])
                        
                        push!(all_distances, dist)
                        push!(all_me_values, me_val)
                    end
                end
                
                # Plot the data
                plot_ME_density_heatmap(
                    all_distances, all_me_values, all_mutations; 
                    alpha=alpha, 
                    y_limit=y_limit
                )
            else
                println("No saved data found at $data_file, need to calculate first.")
                println("Running full analysis instead...")
                # Fall back to full analysis
                plot_only = false
            end
        end
        
        if !plot_only
            # Calculate and create combined plot with all mutations
            all_distances, all_me_values = calculate_combined_ME_data(
                all_mutations, 
                data_path, 
                alpha; 
                num_rounds=10,  # Use fewer rounds for faster processing with many mutations
                verbose=false,  # Less verbose output
                save_data=true
            )
            
            plot_ME_density_heatmap(
                all_distances, all_me_values, all_mutations; 
                alpha=alpha, 
                y_limit=y_limit
            )
        end
    else
        println("No mutations found in directory: $data_path")
    end
end

# Example functions to use:

# 1. Calculate ME data and save it to a file, then plot
# run_full_dataset_analysis("data", 1.0)

# 2. Only plot from previously saved data (much faster)
run_full_dataset_analysis("data", 3.5, plot_only=true)  # y_limit now auto-adapts to data

# 3. Run with a small set of mutations
# run_combined_ME_analysis()

# 4. Separate calculation and plotting
# function example_separate_calc_and_plot()
#     data_path = "data"
#     alpha = 1.0
    
#     # First calculate and save the data
#     all_distances, all_me_values = calculate_combined_ME_data(
#         ["a140e", "a140g", "a140v"], 
#         data_path, 
#         alpha
#     )
    
#     # Then plot with custom parameters
#     plot_ME_density_heatmap(all_distances, all_me_values, y_limit=0.15)
# end

# Run with all mutations and save the data
# run_full_dataset_analysis("data", 3.5)
