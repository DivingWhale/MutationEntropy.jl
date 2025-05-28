using CairoMakie

"""
    test_mutants_with_alpha(mutations, data_path, alpha, residue_range, num_rounds, verbose, plot_results)

Calculate Mutation Entropy (ME) for multiple protein mutants.

This function calculates mutation entropy values for a list of protein mutants
based on a given alpha parameter value, and optionally plots the results.

# Arguments
- `mutations::Vector{String}`: List of mutation identifiers to analyze
- `data_path::String="data"`: Base directory containing protein data
- `alpha::Float64=2.0`: Exponent parameter for strain calculation
- `residue_range=83:231`: Range of residues to analyze
- `num_rounds::Int=20`: Number of simulation rounds to average over
- `verbose::Bool=true`: Whether to print progress information
- `plot_results::Bool=true`: Whether to generate and save plots
"""
function test_mutants_with_alpha(
    mutations::Vector{String}, 
    data_path::String="data", 
    alpha::Float64=2.0, 
    residue_range=83:231, 
    num_rounds::Int=20, 
    verbose::Bool=true,
    plot_results::Bool=true
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
    
    # Process each mutant
    for (i, mutant) in enumerate(mutations)
        println("Processing mutant $i/$(length(mutations)): $mutant")
        
        # Calculate mutant strain values
        mutant_avg_S_values = MutationEntropy.collect_strains(
            data_path,            # Base data directory
            mutant,               # Mutant name
            alpha=alpha,          # Alpha parameter for strain calculation
            residue_range=residue_range, # Range of residues to analyze
            num_rounds=num_rounds, # Number of rounds to average
            verbose=verbose       # Show progress information
        )
        
        # Calculate Mutation Entropy
        mutant_ME = MutationEntropy.calculate_ME(mutant_avg_S_values, wt_avg_S_values)
        
        # Plot results if requested
        if plot_results
            println("Creating plot for $mutant...")
            fig = MutationEntropy.plot_MEvsDist(mutant_ME, data_path, mutant)
            
            # Create figs directory if it doesn't exist
            mkpath("figs")
            
            # Save plot
            save_path = joinpath("figs", "ME_$(mutant)_alpha_$(alpha).png")
            save(save_path, fig)
            println("Plot saved to: $save_path")
        end
    end
    
    println("Completed testing all mutants with alpha = $alpha")
end

# Example: Test multiple mutants
# Define the list of mutations to test
mutations = ["a140e", "a140g", "a140v", "a142c", "a142f", "a142g", "a142v", "a151g", "a151t", "a151v"]
# Set alpha value
alpha = 1.0
# Run the test
test_mutants_with_alpha(mutations, "data", alpha)
