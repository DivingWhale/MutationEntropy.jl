"""
    calculate_mse_per_position(merged_df::DataFrame, ddG_corrected::Vector{Float64})

Calculate Mean Squared Error (MSE) for each position (residue site) based on corrected ddG values.

Returns a DataFrame with columns: position, mse, n_mutations, and a result DataFrame with position, error, corrected_ddG
"""
function calculate_mse_per_position(merged_df::DataFrame, ddG_corrected::Vector{Float64})
    # Add position column if it doesn't exist
    if !(:position in names(merged_df))
        merged_df[!, :position] = [parse_mutation_position(m) for m in merged_df.mutant]
    end

    # Calculate error (difference)
    errors = ddG_corrected .- merged_df.ddG

    # Add error column to DataFrame
    result_df = copy(merged_df)
    result_df[!, :error] = errors
    result_df[!, :corrected_ddG] = ddG_corrected

    # Calculate MSE per position
    mse_df = combine(DataFrames.groupby(result_df, :position),
        :error => (x -> mean(x.^2)) => :mse,
        :error => length => :n_mutations
    )

    return mse_df, result_df
end

"""
    plot_swarm_by_position(result_df::DataFrame, predictor_name::String, output_path::String, A::Float64)

Create swarm plots (beeswarm) showing error distribution for each mutated position.
Handles >20 sites by creating multiple subplots with 20 sites per row.

# Arguments
- `result_df`: DataFrame containing position and error columns
- `predictor_name`: Name of the predictor (e.g., "rosetta", "foldx")
- `output_path`: Path to save the plot
- `A`: The A parameter value used for correction
"""
function plot_swarm_by_position(result_df::DataFrame, predictor_name::String, output_path::String, A::Float64)
    # Check for empty data
    if isempty(result_df)
        @warn "Empty result DataFrame, skipping plot: $output_path"
        return
    end
    
    mkpath(dirname(output_path))

    # Get all positions and sort them
    all_positions = sort(unique(result_df.position))
    n_positions = length(all_positions)
    
    if n_positions == 0
        @warn "No positions found in result DataFrame, skipping plot: $output_path"
        return
    end
    
    sites_per_row = 20
    n_rows = ceil(Int, n_positions / sites_per_row)

    # Calculate figure size - much larger to accommodate bigger fonts
    subplot_width = 3200  # Much larger width for 20 sites with big fonts
    subplot_height = 1200  # Much larger height per row
    fig = Figure(size=(subplot_width, subplot_height * n_rows))

    # Calculate global y-axis range across all positions for consistent scaling
    all_errors = result_df.error
    y_min = minimum(all_errors)
    y_max = maximum(all_errors)
    y_range = y_max - y_min
    y_padding = y_range * 0.1  # 10% padding
    y_limits = (y_min - y_padding, y_max + y_padding)

    # Create subplots
    for row_idx in 1:n_rows
        start_idx = (row_idx - 1) * sites_per_row + 1
        end_idx = min(row_idx * sites_per_row, n_positions)
        positions_subset = all_positions[start_idx:end_idx]

        # Create subplot with consistent y-axis limits and larger fonts
        ax = Axis(fig[row_idx, 1],
            title="Sites $(positions_subset[1])-$(positions_subset[end])",
            xlabel="Residue Position",
            ylabel="ddG_corrected - ddG_exp",
            limits=(nothing, y_limits),  # Use global y-limits for consistency
            titlesize=28,         # Much larger title font
            xlabelsize=24,        # Larger x-axis label font
            ylabelsize=24,        # Larger y-axis label font
            xticklabelsize=20,    # Larger x-axis tick labels
            yticklabelsize=20)    # Larger y-axis tick labels

        # Set y-axis to clip violins that extend beyond limits
        ax.ytickformat = xs -> [string(round(x, digits=1)) for x in xs]

        # Filter data for this subset of positions
        subset_df = result_df[result_df.position .∈ Ref(positions_subset), :]

        # Create swarm plot using SwarmMakie
        beeswarm!(ax, subset_df.position, subset_df.error,
            color=RGBf(0.0431, 0.6471, 0.8745),  # #0BA6DF in RGB
            markersize=10,
            strokecolor=:black,
            strokewidth=1.0)

        # Add zero line for reference - thicker line
        hlines!(ax, 0.0, color=:black, linestyle=:dash, linewidth=2.0)

        # Ensure violins stay within the specified y-limits
        # This must be done after creating all plot elements
        ylims!(ax, y_limits[1], y_limits[2])

        # Rotate x-axis labels for better readability
        ax.xticklabelrotation = 45
    end

    # Adjust layout with more spacing for larger fonts
    rowgap!(fig.layout, 40)   # More space between rows
    colgap!(fig.layout, 30)   # More space between columns

    save(output_path, fig)
    println("Swarm plot saved to: $output_path")
end

function perform_correlation_analysis(merged_df::DataFrame, A_values::Vector{Float64}, predictor_name::String, output_dir::String, prefix::String, target_alpha::Float64, target_rho::Float64, if_normalize::Bool, nearby_normalize::Bool; generate_plots::Bool=false)
    # --- Linear Scaling of Predictor ddG ---
    predictor_col = Symbol(predictor_name * "_ddG")
    scaled_predictor_col = Symbol(predictor_name * "_ddG_scaled")
    scaling_model = lm(Term(:ddG) ~ Term(predictor_col), merged_df)
    intercept = coef(scaling_model)[1]
    slope = coef(scaling_model)[2]
    merged_df[!, scaled_predictor_col] = slope .* merged_df[!, predictor_col] .+ intercept

    # --- Calculate Correlations ---
    # Correlation between Predictor and Experimental ddG (the baseline)
    cor_predictor_exp = cor(merged_df[!, scaled_predictor_col], merged_df.ddG)

    # Correlation for self-mutation ddS vs. Experimental ddG
    cor_mutant_self = cor(merged_df.mutant_ddS, merged_df.ddG)

    # Correlation for nearby residues ddS vs. Experimental ddG
    cor_nearby = cor(merged_df.nearby_ddS, merged_df.ddG)

    # Calculate the difference between experimental and Predictor ddG
    merged_df.ddG_diff = merged_df.ddG .- merged_df[!, scaled_predictor_col]

    # Correlation for self-mutation ddS vs. the ddG difference
    cor_mutant_self_diff = cor(merged_df.mutant_ddS, merged_df.ddG_diff)

    # Correlation for nearby residues ddS vs. the ddG difference
    cor_nearby_diff = cor(merged_df.nearby_ddS, merged_df.ddG_diff)

    # --- Save correlations to file ---
    mkpath(output_dir)
    
    output_file = joinpath(output_dir, "$(prefix)$(predictor_name)_correlations_alpha_$(target_alpha)_rho_$(target_rho).txt")
    open(output_file, "w") do io
        println(io, "$(titlecase(predictor_name)) Correlation Analysis Results")
        println(io, "Alpha: $target_alpha, Rho: $target_rho")
        println(io, "Normalize: $if_normalize, Nearby_normalize: $nearby_normalize")
        println(io, "Number of mutants: $(nrow(merged_df))")
        println(io, "=" ^ 50)
        println(io, "$(titlecase(predictor_name)) vs Experimental ddG: $(round(cor_predictor_exp, digits=4))")
        println(io, "Self ddS vs Experimental ddG: $(round(cor_mutant_self, digits=4))")
        println(io, "Nearby ddS vs Experimental ddG: $(round(cor_nearby, digits=4))")
        println(io, "Self ddS vs ddG_diff: $(round(cor_mutant_self_diff, digits=4))")
        println(io, "Nearby ddS vs ddG_diff: $(round(cor_nearby_diff, digits=4))")
        println(io, "")
        println(io, "Combined correlations ($(titlecase(predictor_name)) + A*ddS):")
        
        plot_dir = joinpath(output_dir, "scatter_plots", "rho_$(target_rho)")
        if generate_plots
            mkpath(plot_dir)
        end

        for A in A_values
            # For self ddS
            combined_self = merged_df[!, scaled_predictor_col] .+ A .* merged_df.mutant_ddS
            cor_self = cor(merged_df.ddG, combined_self)
            
            # Calculate contribution percentages for self ddS
            predictor_contrib_self = abs.(merged_df[!, scaled_predictor_col])
            entropy_contrib_self = abs.(A .* merged_df.mutant_ddS)
            total_contrib_self = predictor_contrib_self .+ entropy_contrib_self
            predictor_pct_self = mean(predictor_contrib_self ./ total_contrib_self) * 100
            entropy_pct_self = mean(entropy_contrib_self ./ total_contrib_self) * 100

            # For nearby ddS
            combined_nearby = merged_df[!, scaled_predictor_col] .+ A .* merged_df.nearby_ddS
            cor_nearby_A = cor(merged_df.ddG, combined_nearby)
            
            # Calculate contribution percentages for nearby ddS
            predictor_contrib_nearby = abs.(merged_df[!, scaled_predictor_col])
            entropy_contrib_nearby = abs.(A .* merged_df.nearby_ddS)
            total_contrib_nearby = predictor_contrib_nearby .+ entropy_contrib_nearby
            predictor_pct_nearby = mean(predictor_contrib_nearby ./ total_contrib_nearby) * 100
            entropy_pct_nearby = mean(entropy_contrib_nearby ./ total_contrib_nearby) * 100
            
            println(io, "A=$A: Self=$(round(cor_self, digits=4)) [$(titlecase(predictor_name)):$(round(predictor_pct_self,digits=1))%, Entropy:$(round(entropy_pct_self,digits=1))%], Nearby=$(round(cor_nearby_A, digits=4)) [$(titlecase(predictor_name)):$(round(predictor_pct_nearby,digits=1))%, Entropy:$(round(entropy_pct_nearby,digits=1))%]")

            if generate_plots
                plot_output_path = joinpath(plot_dir, "$(prefix)$(predictor_name)_scatter_A_$(A).png")
                plot_title = "$(titlecase(predictor_name)) ddG Correction"
                plot_correction_scatter(merged_df.ddG, merged_df[!, scaled_predictor_col], combined_self, cor_self, A, plot_output_path, plot_title)
            end
        end
    end
    
    println("Correlations saved to: $output_file")
end