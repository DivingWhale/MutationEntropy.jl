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

function plot_simple_correlation(x_data::AbstractVector, y_data::AbstractVector, x_label::String, y_label::String, title_text::String, output_path::String)
    if length(x_data) != length(y_data)
        error("Input vectors must have the same length.")
    end
    
    if isempty(x_data)
        @warn "Input data is empty, skipping plot generation."
        return
    end

    # Calculate Pearson correlation
    pcc = cor(x_data, y_data)
    
    # Perform linear regression
    lin_reg = lm(@formula(y ~ x), DataFrame(x=x_data, y=y_data))
    y_fit = predict(lin_reg)

    fig = Figure(size=(800, 600))
    ax = Axis(fig[1, 1], 
              xlabel=x_label, 
              ylabel=y_label, 
              title=title_text,
              titlesize=20)

    scatter!(ax, x_data, y_data, label="Data points")
    lines!(ax, x_data, y_fit, color=:red, label="Regression line")
    
    # Add PCC value to the plot
    text!(ax, "$(round(pcc, digits=3))",
          position=(maximum(x_data), minimum(y_data)),
          align=(:right, :bottom),
          fontsize=18,
          color=:blue)

    axislegend(ax)
    
    # Ensure output directory exists
    mkpath(dirname(output_path))
    
    save(output_path, fig)
    println("Correlation plot saved to: $output_path")
end

function plot_correction_scatter(exp_ddG, pred_ddG, corrected_ddG, pcc, A, output_path, title)
    final_title = "$title (A=$(round(A, digits=3)))"
    
    fig = Figure(size = (800, 600))
    ax = Axis(fig[1, 1], xlabel = "Experimental ddG", ylabel = "Predicted ddG", title = final_title)

    scatter!(ax, exp_ddG, pred_ddG, color = (:grey, 0.7), label = "Predictor ddG")
    scatter!(ax, exp_ddG, corrected_ddG, color = (:orange, 0.7), label = "Corrected ddG")

    min_val = minimum([minimum(exp_ddG), minimum(pred_ddG), minimum(corrected_ddG)])
    max_val = maximum([maximum(exp_ddG), maximum(pred_ddG), maximum(corrected_ddG)])
    lines!(ax, [min_val, max_val], [min_val, max_val], color = :red, linestyle = :dash, label = "y=x")

    axislegend(ax, position = :lt)

    text!(ax, 0.98, 0.05, text = "PCC: $(round(pcc, digits=3))", align = (:right, :bottom), space = :relative, fontsize=16)

    save(output_path, fig)
end