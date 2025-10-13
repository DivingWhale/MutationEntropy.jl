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