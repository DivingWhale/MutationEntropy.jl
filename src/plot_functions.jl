"""
    plot_MEvsDist(ME::Dict{Int, Float64}, datadir::String, mutation::String, round::Int=1)

Plot Mutation Entropy (ME) against distance from the mutation site using cached distance matrix.
Uses CairoMakie for plotting.

# Arguments
- `ME::Dict{Int, Float64}`: Dictionary mapping residue IDs to their mutation entropy values
- `datadir::String`: Base directory containing protein data
- `mutation::String`: Mutation identifier (e.g., "A140E")
- `round::Int=1`: Round number to use for distance matrix (default 1)

# Returns
- The figure object
"""
function plot_MEvsDist(ME::Dict{Int, Float64}, datadir::String, mutation::String, round::Int=1)
    # Load cached distance matrix for the mutation or wild type
    dists = get_dist_map(datadir, mutation, round)
    
    # Extract mutation site residue ID (e.g., get 140 from "A140E")
    mutation_site = parse(Int, match(r"[a-zA-Z](\d+)[a-zA-Z]", mutation).captures[1])
    
    # Extract residue IDs and mutation entropy values
    residues = collect(keys(ME))
    me_values = collect(values(ME))
    
    # Get distances from each residue to the mutation site
    distances = [dists[mutation_site, res] for res in residues]
    
    # Create a figure with CairoMakie
    fig = Figure(size=(800, 600))
    ax = Axis(fig[1, 1], 
              title="Mutation Effect vs Distance from $(mutation) site",
              xlabel="Distance from mutation site (Å)",
              ylabel="Mutation Effect")
    
    # Scatter plot
    scatter!(ax, distances, me_values, 
             color=:orange, 
             markersize=8)
    
    # Add residue IDs as text annotations
    for (i, res) in enumerate(residues)
        text!(ax, distances[i], me_values[i], text="$res", 
              align=(:center, :bottom), 
              offset=(0, 3), 
              fontsize=8)
    end
    
    # Save the plot
    save("figs/ME_vs_Dist_$(mutation)_round$(round).png", fig)
    
    return fig
end

"""
    plot_MEvsDist(ME::Dict{Int, Float64}, dists::Matrix{Float64}, 
                 mutation::String, round::Int=1)

Plot Mutation Entropy (ME) against distance from the mutation site using provided distance matrix.
Uses CairoMakie for plotting.

# Arguments
- `ME::Dict{Int, Float64}`: Dictionary mapping residue IDs to their mutation entropy values
- `dists::Matrix{Float64}`: Distance matrix containing pairwise distances between residues
- `mutation::String`: Mutation identifier (e.g., "A140E")
- `round::Int=1`: Round number (default 1)

# Returns
- The figure object
"""
function plot_MEvsDist(ME::Dict{Int, Float64}, dists::Matrix{Float64}, 
                      mutation::String, round::Int=1)
    # Extract mutation site residue ID (e.g., get 140 from "A140E")
    mutation_site = parse(Int, match(r"[a-zA-Z](\d+)[a-zA-Z]", mutation).captures[1])
    
    # Extract residue IDs and mutation entropy values
    residues = collect(keys(ME))
    me_values = collect(values(ME))
    
    # Get distances from each residue to the mutation site
    distances = [dists[mutation_site, res] for res in residues]
    
    # Create a figure with CairoMakie
    fig = Figure(resolution=(800, 600))
    ax = Axis(fig[1, 1], 
              title="Mutation Effect vs Distance from $(mutation) site",
              xlabel="Distance from mutation site (Å)",
              ylabel="Mutation Effect")
    
    # Scatter plot
    scatter!(ax, distances, me_values, 
             color=:orange, 
             markersize=8)
    
    # Add residue IDs as text annotations
    for (i, res) in enumerate(residues)
        text!(ax, distances[i], me_values[i], text="$res", 
              align=(:center, :bottom), 
              offset=(0, 3), 
              fontsize=8)
    end
    
    # Save the plot
    save("ME_vs_Dist_$(mutation)_round$(round).png", fig)
    
    return fig
end
