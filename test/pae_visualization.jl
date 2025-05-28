using MutationEntropy
using CairoMakie
using Statistics
using LinearAlgebra  # 用于向量计算
using Glob  # 用于文件查找

"""
    plot_pae_matrix(pae::Matrix{Float64}, mutation::String, round_num::Int=1)

Plot the PAE (Predicted Aligned Error) matrix as a heatmap for a specific mutation.

# Arguments
- `pae::Matrix{Float64}`: PAE matrix data
- `mutation::String`: Mutation identifier (e.g., "l120e")
- `round_num::Int=1`: Round number to display in the title

# Returns
- The figure object
"""
function plot_pae_matrix(pae::Matrix{Float64}, mutation::String, round_num::Int=1)
    # Create a new figure
    fig = Figure(size=(800, 700))
    
    # Create an axis with appropriate labels
    ax = Axis(fig[1, 1],
              title="PAE Matrix for $(mutation) (Round $round_num)",
              xlabel="Residue Index",
              ylabel="Residue Index")
    
    # Create a heatmap of the PAE matrix
    hm = heatmap!(ax, pae,
                  colormap=:viridis,  # Reversed viridis (darker = higher error)
                  colorrange=(0, 30))   # Range for the color scale (0-30 Å is typical)
    
    # Add a colorbar
    cbar = Colorbar(fig[1, 2], hm, label="Predicted Aligned Error (Å)")
    
    return fig
end

"""
    plot_pae_distribution(pae::Matrix{Float64}, mutation::String, round_num::Int=1)

Plot a histogram of the PAE (Predicted Aligned Error) values for a specific mutation.

# Arguments
- `pae::Matrix{Float64}`: PAE matrix data
- `mutation::String`: Mutation identifier (e.g., "l120e")
- `round_num::Int=1`: Round number to display in the title

# Returns
- The figure object
"""
function plot_pae_distribution(pae::Matrix{Float64}, mutation::String, round_num::Int=1)
    # Flatten the matrix to get all PAE values
    pae_values = vec(pae)
    
    # Filter out extremely large values for better visualization
    max_pae = min(maximum(pae_values), 30.0)
    pae_values = filter(x -> x <= max_pae, pae_values)
    
    # Create a new figure
    fig = Figure(size=(800, 600))
    
    # Create an axis with appropriate labels
    ax = Axis(fig[1, 1],
              title="PAE Distribution for $(mutation) (Round $round_num)",
              xlabel="Predicted Aligned Error (Å)",
              ylabel="Frequency")
    
    # Create a histogram
    hist!(ax, pae_values, bins=50, color=:skyblue, strokewidth=1, strokecolor=:black)
    
    # Add statistics as text
    median_val = median(pae_values)
    mean_val = mean(pae_values)
    
    text!(ax, 0.7, 0.9, text="Median: $(round(median_val, digits=2)) Å\nMean: $(round(mean_val, digits=2)) Å",
          align=(:left, :top), space=:relative, color=:black)
    
    return fig
end

"""
    compare_pae_matrices(pae1::Matrix{Float64}, label1::String, 
                         pae2::Matrix{Float64}, label2::String)

Compare two PAE matrices side by side with a difference map.

# Arguments
- `pae1::Matrix{Float64}`: First PAE matrix
- `label1::String`: Label for the first matrix (e.g., "Wild Type")
- `pae2::Matrix{Float64}`: Second PAE matrix
- `label2::String`: Label for the second matrix (e.g., "Mutant")

# Returns
- The figure object
"""
function compare_pae_matrices(pae1::Matrix{Float64}, label1::String, 
                             pae2::Matrix{Float64}, label2::String)
    # Check if matrices have the same dimensions
    if size(pae1) != size(pae2)
        error("PAE matrices must have the same dimensions for comparison")
    end
    
    # Calculate the difference matrix (mutant - wildtype)
    diff_matrix = pae2 .- pae1
    
    # Create a new figure with 3 subplots
    fig = Figure(size=(1200, 400))
    
    # Create three axes
    ax1 = Axis(fig[1, 1], title=label1, xlabel="Residue Index", ylabel="Residue Index")
    ax2 = Axis(fig[1, 2], title=label2, xlabel="Residue Index")
    ax3 = Axis(fig[1, 3], title="Difference ($(label2) - $(label1))", xlabel="Residue Index")
    
    # Create heatmaps
    hm1 = heatmap!(ax1, pae1, colormap=:viridis, colorrange=(0, 30))
    hm2 = heatmap!(ax2, pae2, colormap=:viridis, colorrange=(0, 30))
    
    # For the difference, use a diverging colormap centered at 0
    max_diff = max(abs(minimum(diff_matrix)), abs(maximum(diff_matrix)))
    hm3 = heatmap!(ax3, diff_matrix, colormap=:balance, colorrange=(-max_diff, max_diff))
    
    # Add colorbars
    cbar1 = Colorbar(fig[1, 4], hm1, label="PAE (Å)")
    cbar2 = Colorbar(fig[1, 5], hm3, label="ΔPAE (Å)")
    
    return fig
end

"""
    plot_pae_vs_distance(pae::Matrix{Float64}, mutation::String, data_path::String, round_num::Int=1; highlight_cutoff::Float64=13.0)

Plot the PAE (Predicted Aligned Error) values against distance from the mutation site.

# Arguments
- `pae::Matrix{Float64}`: PAE matrix data
- `mutation::String`: Mutation identifier (e.g., "l120e")
- `data_path::String`: Path to data directory containing PDB files
- `round_num::Int=1`: Round number to display in the title
- `highlight_cutoff::Float64=13.0`: Distance cutoff to highlight nearby residues (in Å)

# Returns
- The figure object
"""
function plot_pae_vs_distance(pae::Matrix{Float64}, mutation::String, data_path::String, round_num::Int=1; highlight_cutoff::Float64=13.0)
    # Extract mutation site number
    mut_site = parse(Int, match(r"[a-zA-Z](\d+)[a-zA-Z]", mutation).captures[1])
    
    # Get the PAE values for the mutation site row (interactions from mutation site to all other residues)
    mut_site_row = pae[mut_site, :]
    
    # Get the distance matrix using MutationEntropy's function
    println("Calculating distance matrix...")
    d_matrix = MutationEntropy.get_dist_map(data_path, mutation, round_num)
    
    # Calculate distances from mutation site to all other residues
    distances = Float64[]
    pae_values = Float64[]
    residue_indices = Int[]
    
    # Process all residues
    for i in 1:length(mut_site_row)
        if i != mut_site && i <= size(d_matrix, 2)
            dist = d_matrix[mut_site, i]
            push!(distances, dist)
            push!(pae_values, mut_site_row[i])
            push!(residue_indices, i)
        end
    end
    
    # Create figure
    fig = Figure(size=(800, 600))
    
    # Use MutationEntropy's function to find residues within the cutoff distance
    nearby_residues = MutationEntropy.find_residues_within_distance(mut_site, d_matrix, distance=highlight_cutoff)
    
    # Create figure
    fig = Figure(size=(900, 600))
    
    # Create scatter plot
    ax = Axis(fig[1, 1],
              title="PAE vs Distance from Mutation Site $(mutation) (Round $(round_num))",
              xlabel="Distance from residue $(mut_site) (Å)",
              ylabel="Predicted Aligned Error (Å)")
    
    # Mark points within the cutoff distance with a different color
    nearby_indices = findall(i -> i ∈ nearby_residues, residue_indices)
    distant_indices = findall(i -> !(i ∈ nearby_residues), residue_indices)
    
    # Create scatter plots with different colors for nearby and distant residues
    if !isempty(distant_indices)
        scatter!(ax, distances[distant_indices], pae_values[distant_indices], 
                markersize=8,
                color=:lightgray,
                alpha=0.6,
                label="d > $(highlight_cutoff)Å")
    end
    
    if !isempty(nearby_indices)
        scatter!(ax, distances[nearby_indices], pae_values[nearby_indices], 
                markersize=8,
                color=:darkblue,
                alpha=0.8,
                label="d ≤ $(highlight_cutoff)Å")
    end
    
    # Add trend line (smoothed average)
    sorted_indices = sortperm(distances)
    sorted_distances = distances[sorted_indices]
    sorted_pae = pae_values[sorted_indices]
    
    # Use a simple moving average for the trend line
    window_size = max(5, div(length(sorted_distances), 10))
    trend_distances = Float64[]
    trend_pae = Float64[]
    
    for i in window_size:length(sorted_distances)-window_size
        push!(trend_distances, sorted_distances[i])
        push!(trend_pae, mean(sorted_pae[i-window_size+1:i+window_size-1]))
    end
    
    lines!(ax, trend_distances, trend_pae, 
          color=:red, 
          linewidth=3, 
          label="Trend")
    
    # Add vertical line at the cutoff distance
    vlines!(ax, highlight_cutoff, 
           color=:red, 
           linestyle=:dash, 
           label="Cutoff ($(highlight_cutoff)Å)")
    
    # Add legend
    Legend(fig[1, 2], ax, "Legend", framevisible=true)
    
    # Add annotations for some notable points
    # Find highest PAE values
    top_indices = sortperm(pae_values, rev=true)[1:min(5, length(pae_values))]
    
    for idx in top_indices
        text!(ax, distances[idx], pae_values[idx] + 0.5, 
              text="$(residue_indices[idx])",
              fontsize=12,
              align=(:center, :bottom))
    end
    
    # Calculate correlations and statistics
    overall_correlation = cor(distances, pae_values)
    
    # Statistics for nearby residues
    nearby_stats = ""
    if !isempty(nearby_indices)
        nearby_correlation = cor(distances[nearby_indices], pae_values[nearby_indices])
        nearby_mean = mean(pae_values[nearby_indices])
        nearby_stats = "Within $(highlight_cutoff)Å:\n" *
                       "  Mean PAE: $(round(nearby_mean, digits=2)) Å\n" *
                       "  Correlation: $(round(nearby_correlation, digits=2))\n"
    end
    
    # Add statistical information
    text!(ax, 0.05, 0.95, 
          text="Overall Statistics:\n" *
               "  Correlation: $(round(overall_correlation, digits=2))\n" *
               "  Mean PAE: $(round(mean(pae_values), digits=2)) Å\n" *
               nearby_stats,
          align=(:left, :top), 
          space=:relative,
          color=:black)
    
    return fig
end

"""
    plot_pae_difference_profile(wt_pae::Matrix{Float64}, mut_pae::Matrix{Float64}, 
                              mutation::String, data_path::String, round_num::Int=1)

Plot the difference in PAE values between wild type and mutant for the mutation site.

# Arguments
- `wt_pae::Matrix{Float64}`: Wild type PAE matrix
- `mut_pae::Matrix{Float64}`: Mutant PAE matrix
- `mutation::String`: Mutation identifier (e.g., "l120e")
- `data_path::String`: Path to data directory containing PDB files
- `round_num::Int=1`: Round number to display in the title

# Returns
- The figure object
"""
function plot_pae_difference_profile(wt_pae::Matrix{Float64}, mut_pae::Matrix{Float64}, 
                                    mutation::String, data_path::String, round_num::Int=1)
    # Extract mutation site number
    mut_site = parse(Int, match(r"[a-zA-Z](\d+)[a-zA-Z]", mutation).captures[1])
    
    # Get residue range
    residues = 1:min(size(wt_pae, 1), size(mut_pae, 1))
    
    # Get the PAE values for the mutation site row
    wt_row = wt_pae[mut_site, :]
    mut_row = mut_pae[mut_site, :]
    
    # Calculate difference (mutant - wild type)
    diff_values = mut_row .- wt_row
    
    # Get distance matrix for highlighting nearby residues
    d_matrix = MutationEntropy.get_dist_map(data_path, mutation, round_num)
    nearby_residues = MutationEntropy.find_residues_within_distance(mut_site, d_matrix, distance=13.0)
    
    # Create a new figure
    fig = Figure(size=(1000, 600))
    
    # Create the top panel for PAE profiles
    ax1 = Axis(fig[1, 1],
              title="PAE Profile for Site $(mut_site) ($(mutation))",
              ylabel="PAE (Å)")
    
    # Create the bottom panel for differences
    ax2 = Axis(fig[2, 1],
              title="PAE Difference (Mutant - Wild Type)",
              xlabel="Residue Number",
              ylabel="ΔPAE (Å)")
    
    # Hide x-axis labels on the top plot
    hidexdecorations!(ax1, grid=false)
    
    # Link x-axes for synchronized zooming/panning
    linkxaxes!(ax1, ax2)
    
    # Plot data in top panel
    lines!(ax1, residues, wt_row, color=:blue, label="Wild Type")
    lines!(ax1, residues, mut_row, color=:red, label="Mutant ($(mutation))")
    
    # Add vertical line at mutation site
    vlines!(ax1, mut_site, color=:black, linestyle=:dash, label="Mutation Site")
    vlines!(ax2, mut_site, color=:black, linestyle=:dash)
    
    # Plot data in bottom panel
    stem!(ax2, residues, diff_values)
    
    # Add horizontal line at zero for difference plot
    hlines!(ax2, 0, color=:gray, linestyle=:solid)
    
    # Highlight nearby residues 
    for res in nearby_residues
        if res != mut_site && res <= length(residues)
            # Draw vertical spans to highlight these residues
            vspan!(ax1, res-0.5, res+0.5, color=(:yellow, 0.2))
            vspan!(ax2, res-0.5, res+0.5, color=(:yellow, 0.2))
        end
    end
    
    # Add legend
    Legend(fig[1, 2], ax1)
    
    # Add annotations about significant changes
    # Find top 3 increases and decreases
    num_to_highlight = 3
    sorted_diffs = sort(collect(zip(residues, diff_values)), by=x->x[2])
    
    top_decreases = sorted_diffs[1:min(num_to_highlight, length(sorted_diffs))]
    top_increases = reverse(sorted_diffs[max(1, length(sorted_diffs)-num_to_highlight+1):end])
    
    # Add text boxes with statistics
    textbox = "Mutation: $(mutation)\n" *
              "Site: $(mut_site)\n\n" *
              "Largest PAE Increases:\n"
    
    for (res, diff) in top_increases
        textbox *= "  Res $(res): +$(round(diff, digits=2)) Å\n"
    end
    
    textbox *= "\nLargest PAE Decreases:\n"
    for (res, diff) in top_decreases
        textbox *= "  Res $(res): $(round(diff, digits=2)) Å\n"
    end
    
    # Calculate mean changes
    mean_change_all = mean(diff_values)
    mean_change_nearby = mean([diff_values[i] for i in nearby_residues if i <= length(diff_values)])
    
    textbox *= "\nMean Changes:\n" *
               "  All residues: $(round(mean_change_all, digits=2)) Å\n" *
               "  Nearby (<13Å): $(round(mean_change_nearby, digits=2)) Å"
    
    # Add text annotation
    text!(ax1, 0.98, 0.98, text=textbox,
          align=(:right, :top), 
          space=:relative,
          color=:black,
          fontsize=12)
    
    return fig
end

# Main execution code
data_path = "data"
mutant = "a140g"
round_num = 1

# Create directory for figures if it doesn't exist
if !isdir("figs")
    mkdir("figs")
end

# Read PAE matrices
println("Reading PAE matrix for mutant $(mutant)...")
mutant_pae = MutationEntropy.read_pae(data_path, mutant, round_num)

println("Reading PAE matrix for wild type...")
wt_pae = MutationEntropy.read_pae(data_path, "thermonulease", round_num)

# 1. Plot the mutant PAE matrix
println("Generating PAE matrix heatmap...")
pae_fig = plot_pae_matrix(mutant_pae, mutant, round_num)
save(joinpath("figs", "$(mutant)_pae_matrix.png"), pae_fig)
println("Saved PAE matrix heatmap to figs/$(mutant)_pae_matrix.png")

# 2. Plot the PAE value distribution
println("Generating PAE distribution histogram...")
dist_fig = plot_pae_distribution(mutant_pae, mutant, round_num)
save(joinpath("figs", "$(mutant)_pae_distribution.png"), dist_fig)
println("Saved PAE distribution to figs/$(mutant)_pae_distribution.png")

# 3. Compare with wild type PAE
println("Generating PAE comparison visualization...")
comp_fig = compare_pae_matrices(wt_pae, "Wild Type", mutant_pae, mutant)
save(joinpath("figs", "$(mutant)_vs_wt_pae.png"), comp_fig)
println("Saved comparison visualization to figs/$(mutant)_vs_wt_pae.png")

# 4. Plot PAE vs distance from mutation site
println("Generating PAE vs distance visualization...")
dist_pae_fig = plot_pae_vs_distance(mutant_pae, mutant, data_path, round_num, highlight_cutoff=13.0)
save(joinpath("figs", "$(mutant)_pae_vs_distance.png"), dist_pae_fig)
println("Saved PAE vs distance visualization to figs/$(mutant)_pae_vs_distance.png")

# 5. Plot PAE difference profile for the mutation site
println("Generating PAE difference profile...")
diff_profile_fig = plot_pae_difference_profile(wt_pae, mutant_pae, mutant, data_path, round_num)
save(joinpath("figs", "$(mutant)_pae_difference_profile.png"), diff_profile_fig)
println("Saved PAE difference profile to figs/$(mutant)_pae_difference_profile.png")

# 6. Optionally extract data about the mutation site
mut_site = parse(Int, match(r"[a-zA-Z](\d+)[a-zA-Z]", mutant).captures[1])
println("\nMutation site: $mut_site")

# Extract PAE values for the mutation site row/column
mut_site_row = mutant_pae[mut_site, :]
mut_site_col = mutant_pae[:, mut_site]

# Get mean PAE for the mutation site (excluding self-interaction)
mean_pae = mean([mut_site_row[i] for i in 1:length(mut_site_row) if i != mut_site])
println("Mean PAE from residue $mut_site to all others: $(round(mean_pae, digits=2)) Å")

# Calculate difference from wild type specifically at mutation site
wt_mut_site_row = wt_pae[mut_site, :]
wt_mut_site_col = wt_pae[:, mut_site]

println("PAE difference at mutation site $mut_site comparing to wild type:")
delta_mean_pae = mean([mut_site_row[i] - wt_mut_site_row[i] for i in 1:length(mut_site_row) if i != mut_site])
println("Mean ΔPAE: $(round(delta_mean_pae, digits=2)) Å")

# If the result is positive, the mutant has higher error (less confident)
if delta_mean_pae > 0
    println("The mutation increases structural uncertainty around position $mut_site.")
else
    println("The mutation decreases structural uncertainty around position $mut_site.")
end
