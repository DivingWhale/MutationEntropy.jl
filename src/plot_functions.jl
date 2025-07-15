function plot_MEvsDist(ME::Dict{Int, Float64}, mutation::String, dists::Matrix{Float64}; datadir::String="data", round::Int=1)
    mutation_site = parse(Int, match(r"[a-zA-Z](\d+)[a-zA-Z]", mutation).captures[1])
    low_plddt_residues = get_low_plddt_residues(mutation, round, datadir)
    residues = collect(keys(ME))
    residues = filter(r -> !(r in low_plddt_residues), residues)
    me_values = [ME[r] for r in residues]
    distances = [dists[mutation_site, res] for res in residues]
    
    fig = Figure(size=(800, 600))
    ax = Axis(fig[1, 1], 
              title="Mutation Effect vs Distance from $(mutation) site",
              xlabel="Distance from mutation site (Å)",
              ylabel="Mutation Effect")
    
    scatter!(ax, distances, me_values, 
             color=:orange, 
             markersize=8)
    
    for (i, res) in enumerate(residues)
        text!(ax, distances[i], me_values[i], text="$res", 
              align=(:center, :bottom), 
              offset=(0, 3), 
              fontsize=8)
    end
    
    return fig
end

function plot_MEvsDist(ME::Dict{Int, Float64}, datadir::String, mutation::String, round::Int=1)
    dists = get_dist_map(datadir, mutation, round)
    return plot_MEvsDist(ME, mutation, dists; datadir=datadir, round=round)
end