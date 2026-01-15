"""
PCC Plotting Functions

This module provides functions for generating PCC vs A plots.
"""

using CairoMakie
using Colors

export plot_pcc_vs_A, get_baseline_correlation

"""
    plot_pcc_vs_A(; kwargs...)

Plot correlation coefficients (PCC) vs A values for different parameter combinations.
Creates one plot per rho value, with lines representing different normalization settings.

# Keyword Arguments
- `methods::Vector{String}`: Predictors to plot (default: all four)
- `alpha::Float64`: Alpha value (default: 0.0)
- `rho_values`: Rho values to iterate over
- `threshold::String`: Threshold value (default: "70")
- `base_dir::String`: Base directory for figures (default: "figs")
- `ylims`: Y-axis limits tuple or nothing for auto
- `title_fontsize`, `label_fontsize`, `ticklabel_fontsize`, `legend_fontsize`: Font sizes
- `legend_position`: Legend position symbol (e.g., :rb, :rt)
"""
function plot_pcc_vs_A(;
    methods::Vector{String} = ["foldx", "pythia", "rosetta", "stabilityoracle"],
    alpha::Float64 = 0.0,
    rho_values = [8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0],
    threshold::String = "70",
    base_dir::String = "figs",
    ylims = nothing,
    title_fontsize::Union{Int, Nothing} = nothing,
    label_fontsize::Union{Int, Nothing} = nothing,
    ticklabel_fontsize::Union{Int, Nothing} = nothing,
    legend_fontsize::Union{Int, Nothing} = nothing,
    legend_position::Union{Symbol, Nothing} = nothing,
)
    
    for method in methods
        println("Processing $method...")
        
        # Set up output directory based on method
        base_output_dir = if method == "foldx"
            joinpath(base_dir, "entropy_foldx", threshold)
        elseif method == "pythia"
            joinpath(base_dir, "entropy_pythia", threshold)
        elseif method == "rosetta"
            joinpath(base_dir, "entropy", threshold)
        elseif method == "stabilityoracle"
            joinpath(base_dir, "entropy_stabilityoracle", threshold)
        else
            joinpath(base_dir, "entropy_$method", threshold)
        end
        
        for rho in rho_values
            try
                _plot_single_pcc_vs_A(
                    method, base_output_dir, alpha, rho,
                    ylims, title_fontsize, label_fontsize, 
                    ticklabel_fontsize, legend_fontsize, legend_position
                )
            catch e
                println("  ✗ Error processing rho=$rho: $e")
            end
        end
    end
end

"""
Internal function to plot a single PCC vs A figure.
"""
function _plot_single_pcc_vs_A(
    method::String, 
    base_output_dir::String, 
    alpha::Float64, 
    rho::Float64,
    ylims,
    title_fontsize,
    label_fontsize,
    ticklabel_fontsize,
    legend_fontsize,
    legend_position
)
    # Read correlation data for all three normalization combinations
    data_combinations = [
        (false, false, ""),       # if_normalize=false, nearby_normalize=false
        (true, false, "n_"),      # if_normalize=true, nearby_normalize=false
        (true, true, "n_n2_")     # if_normalize=true, nearby_normalize=true
    ]
    
    fig = Figure(size = (1400, 900), fontsize = 16)
    ax = Axis(fig[1, 1],
        xlabel = "A value",
        ylabel = "Correlation Coefficient (PCC)",
        title = "$(uppercase(method)) PCC vs A - Alpha=$(alpha), Rho=$(rho)",
        limits = ylims === nothing ? (nothing, nothing, 0, nothing) : (nothing, nothing, ylims[1], ylims[2]),
        titlesize = title_fontsize === nothing ? 16 : title_fontsize,
        xlabelsize = label_fontsize === nothing ? 16 : label_fontsize,
        ylabelsize = label_fontsize === nothing ? 16 : label_fontsize,
        xticklabelsize = ticklabel_fontsize === nothing ? 16 : ticklabel_fontsize,
        yticklabelsize = ticklabel_fontsize === nothing ? 16 : ticklabel_fontsize
    )
    
    colors = [colorant"#03A6A1", colorant"#FF7D29", colorant"#932F67"]
    labels = ["No normalization", "ddS normalized", "Both normalized"]
    linestyles = [:solid, :dash]  # solid for self, dash for nearby
    
    for (i, (if_normalize, nearby_normalize, prefix)) in enumerate(data_combinations)
        correlation_file = joinpath(base_output_dir, 
            "$(prefix)$(method)_correlations_alpha_$(alpha)_rho_$(rho).txt")
        
        if isfile(correlation_file)
            # Plot both self and nearby correlations
            for (j, corr_type) in enumerate(["self", "nearby"])
                if i == 3 && corr_type == "self"
                    continue # Skip "Both normalized (self)"
                end
                A_values, pcc_values = read_correlation_data(correlation_file, corr_type)
                
                if !isempty(A_values)
                    label_text = "$(labels[i]) ($(corr_type))"
                    lines!(ax, A_values, pcc_values, 
                        color = colors[i], 
                        linewidth = 2, 
                        linestyle = linestyles[j],
                        label = label_text)
                end
            end
        else
            println("  ⚠ File not found: $correlation_file")
        end
    end
    
    try
        axislegend(ax, position = legend_position === nothing ? :rb : legend_position, 
                   nbanks = 2, labelsize = legend_fontsize === nothing ? 16 : legend_fontsize)
    catch e
        if e isa ArgumentError && occursin("no plots with labels", e.msg)
            println("  No labeled data to create legend for $method rho=$rho")
        else
            rethrow(e)
        end
    end
    
    # Add grid
    ax.xgridvisible = true
    ax.ygridvisible = true
    ax.xgridcolor = (:gray, 0.3)
    ax.ygridcolor = (:gray, 0.3)
    
    # Add horizontal reference lines
    hlines!(ax, [0.0], color = :black, linestyle = :solid, alpha = 0.7, linewidth = 1)
    
    # Add baseline correlation reference line (method vs experimental ddG)
    baseline_pcc = get_baseline_correlation(method, base_output_dir, alpha, rho)
    if baseline_pcc !== nothing
        hlines!(ax, [baseline_pcc], color = :gray, linestyle = :dash, alpha = 0.8, linewidth = 2)
        # Get A range for text positioning
        all_A_values = Float64[]
        for (_, _, prefix) in data_combinations
            correlation_file = joinpath(base_output_dir, 
                "$(prefix)$(method)_correlations_alpha_$(alpha)_rho_$(rho).txt")
            if isfile(correlation_file)
                A_vals, _ = read_correlation_data(correlation_file, "self")
                if !isempty(A_vals)
                    append!(all_A_values, A_vals)
                end
            end
        end
        # Add text label for the baseline
        if !isempty(all_A_values)
            text!(ax, "$(uppercase(method)) baseline: $(round(baseline_pcc, digits=3))", 
                position = (minimum(all_A_values), baseline_pcc + 0.02), 
                align = (:left, :bottom), fontsize = 18, color = :gray)
        end
    end
    
    # Save figure
    mkpath(base_output_dir)
    output_file = joinpath(base_output_dir, "$(method)_pcc_vs_A_combined_rho_$(rho).png")
    save(output_file, fig)
    println("  ✓ Saved: $output_file")
end

"""
    get_baseline_correlation(method, base_output_dir, alpha, rho)

Get baseline correlation (method vs experimental ddG) from correlation file.
"""
function get_baseline_correlation(method::String, base_output_dir::String, alpha::Float64, rho::Float64)
    # Try to find any correlation file for this method/rho combination
    # Priority: n_n2_ > n_ > ""
    prefixes = ["n_n2_", "n_", ""]
    
    for prefix in prefixes
        correlation_file = joinpath(base_output_dir, 
            "$(prefix)$(method)_correlations_alpha_$(alpha)_rho_$(rho).txt")
        
        if isfile(correlation_file)
            try
                lines = readlines(correlation_file)
                for line in lines
                    # Look for the baseline correlation line
                    method_lower = lowercase(method)
                    if occursin(Regex("$(method_lower) vs Experimental ddG:", "i"), line)
                        return parse(Float64, split(line, ":")[2] |> strip)
                    end
                end
            catch e
                continue
            end
        end
    end
    
    return nothing
end
