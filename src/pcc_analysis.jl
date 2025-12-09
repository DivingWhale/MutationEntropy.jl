"""
PCC (Pearson Correlation Coefficient) Analysis Functions

This module provides functions for analyzing correlation results from different predictors,
including best PCC analysis, fixed rho analysis, and result parsing.
"""

using Printf
using Dates

export analyze_pcc_files, analyze_fixed_rho, load_pcc_results, read_correlation_data
export find_best_pcc_parameters, save_pcc_results

"""
    analyze_pcc_files(dir_path; io=stdout)

Analyze PCC results from correlation files in the specified directory.
Supports Rosetta, FoldX, Pythia, and StabilityOracle predictors.

# Arguments
- `dir_path`: Directory containing the correlation result files
- `io`: IO stream for output (default: stdout)

# Returns
- Dict containing analysis results with best parameters for each normalization type
"""
function analyze_pcc_files(dir_path; io=stdout)
    files = filter(f -> endswith(f, ".txt"), readdir(dir_path, join=true))
    
    if isempty(files)
        println(io, "No .txt files found in directory")
        return nothing
    end
    
    # Auto-detect file prefix from first file
    first_file = basename(files[1])
    predictor_prefix = ""
    if occursin("rosetta_correlations", first_file)
        predictor_prefix = "rosetta_correlations"
    elseif occursin("foldx_correlations", first_file)
        predictor_prefix = "foldx_correlations"
    elseif occursin("pythia_correlations", first_file)
        predictor_prefix = "pythia_correlations"
    elseif occursin("stabilityoracle_correlations", first_file)
        predictor_prefix = "stabilityoracle_correlations"
    else
        println(io, "Warning: Unknown file format, cannot detect predictor type")
        return nothing
    end
    
    # Define normalization types based on detected prefix
    normalization_types = [
        (predictor_prefix, "No normalization"),
        ("n_$(predictor_prefix)", "Normalization (n)"),
        ("n_n2_$(predictor_prefix)", "Normalization (n_n2)")
    ]
    
    println(io, "="^80)
    println(io, "Analyzing PCC results in: $dir_path")
    println(io, "Detected predictor: $(replace(predictor_prefix, "_correlations" => ""))")
    println(io, "="^80)
    println(io)
    
    all_results = Dict{String, Any}()
    
    for (prefix, type_name) in normalization_types
        println(io, "Type: $type_name")
        println(io, "-"^80)
        
        # Filter files for this type - use startswith for exact prefix match
        type_files = filter(f -> startswith(basename(f), prefix), files)
        
        # Store results
        results = Dict{Float64, Dict{String, Any}}()
        
        # Store base correlations (will be same for all rho values)
        base_correlations = nothing
        
        for file in type_files
            # Extract rho value
            m = match(r"rho_(\d+\.?\d*)\.txt", basename(file))
            if m === nothing
                continue
            end
            rho = parse(Float64, m.captures[1])
            
            # Read file content
            content = read(file, String)
            lines = split(content, '\n')
            
            # Extract base correlations (only once)
            if base_correlations === nothing
                base_correlations = _extract_base_correlations(lines)
            end
            
            # Extract best PCC for Self and Nearby (excluding A=0)
            max_self = -Inf
            max_nearby = -Inf
            max_self_A = nothing
            max_nearby_A = nothing
            
            for line in lines
                # Match format: "A=xxx: Self=yyy [...], Nearby=zzz [...]"
                m = match(r"A=([-\d.]+):\s+Self=([-\d.]+).*Nearby=([-\d.]+)", line)
                if m !== nothing
                    A = parse(Float64, m.captures[1])
                    # Skip A=0 as it's just the original predictor result
                    if A == 0.0
                        continue
                    end
                    
                    self_val = parse(Float64, m.captures[2])
                    nearby_val = parse(Float64, m.captures[3])
                    
                    if self_val > max_self
                        max_self = self_val
                        max_self_A = A
                    end
                    
                    if nearby_val > max_nearby
                        max_nearby = nearby_val
                        max_nearby_A = A
                    end
                end
            end
            
            results[rho] = Dict(
                "max_self" => max_self,
                "max_self_A" => max_self_A,
                "max_nearby" => max_nearby,
                "max_nearby_A" => max_nearby_A
            )
        end
        
        # Find global best
        best_self_rho = nothing
        best_self_pcc = -Inf
        best_self_A = nothing
        
        best_nearby_rho = nothing
        best_nearby_pcc = -Inf
        best_nearby_A = nothing
        
        for (rho, vals) in sort(collect(results), by=x->x[1])
            if vals["max_self"] > best_self_pcc
                best_self_pcc = vals["max_self"]
                best_self_rho = rho
                best_self_A = vals["max_self_A"]
            end
            
            if vals["max_nearby"] > best_nearby_pcc
                best_nearby_pcc = vals["max_nearby"]
                best_nearby_rho = rho
                best_nearby_A = vals["max_nearby_A"]
            end
        end
        
        # Display base correlations
        if base_correlations !== nothing && !isempty(base_correlations)
            println(io, "\nBase Correlations:")
            # Predictor vs Experimental
            for key in sort(collect(keys(base_correlations)))
                if occursin("_vs_exp", key) && !occursin("ddS", key)
                    predictor = replace(key, "_vs_exp" => "")
                    @printf(io, "  %s vs Experimental ddG: %.4f\n", 
                            predictor, base_correlations[key])
                end
            end
            # ddS correlations
            if haskey(base_correlations, "self_ddS_vs_exp")
                @printf(io, "  Self ddS vs Experimental ddG: %.4f\n", 
                        base_correlations["self_ddS_vs_exp"])
            end
            if haskey(base_correlations, "nearby_ddS_vs_exp")
                @printf(io, "  Nearby ddS vs Experimental ddG: %.4f\n", 
                        base_correlations["nearby_ddS_vs_exp"])
            end
            if haskey(base_correlations, "self_ddS_vs_diff")
                @printf(io, "  Self ddS vs ddG_diff: %.4f\n", 
                        base_correlations["self_ddS_vs_diff"])
            end
            if haskey(base_correlations, "nearby_ddS_vs_diff")
                @printf(io, "  Nearby ddS vs ddG_diff: %.4f\n", 
                        base_correlations["nearby_ddS_vs_diff"])
            end
        end
        
        # Display best results
        println(io, "\nBest Combined Results (Predictor + A*ddS):")
        if best_self_rho !== nothing && best_self_A !== nothing
            @printf(io, "  Best Self:   PCC=%.4f, rho=%.1f, A=%.3f\n", 
                    best_self_pcc, best_self_rho, best_self_A)
        else
            println(io, "  Best Self:   No valid results found")
        end
        
        if best_nearby_rho !== nothing && best_nearby_A !== nothing
            @printf(io, "  Best Nearby: PCC=%.4f, rho=%.1f, A=%.3f\n", 
                    best_nearby_pcc, best_nearby_rho, best_nearby_A)
        else
            println(io, "  Best Nearby: No valid results found")
        end
        
        println(io, "\n" * "="^80)
        println(io)
        
        # Store results
        all_results[type_name] = Dict(
            "best_self_pcc" => best_self_pcc,
            "best_self_rho" => best_self_rho,
            "best_self_A" => best_self_A,
            "best_nearby_pcc" => best_nearby_pcc,
            "best_nearby_rho" => best_nearby_rho,
            "best_nearby_A" => best_nearby_A,
            "base_correlations" => base_correlations
        )
    end
    
    return all_results
end

"""
    analyze_fixed_rho(dir_path, target_rho; io=stdout)

Analyze correlation results for a fixed rho value across all methods and normalization types.

# Arguments
- `dir_path`: Directory containing correlation files
- `target_rho`: The target rho value to analyze
- `io`: IO stream for output (default: stdout)

# Returns
- Dict containing analysis results for the fixed rho
"""
function analyze_fixed_rho(dir_path, target_rho; io=stdout)
    if !isdir(dir_path)
        println(io, "Warning: Directory does not exist: $dir_path")
        return nothing
    end
    
    # Get all txt files
    files = filter(f -> endswith(f, ".txt"), 
                   [joinpath(dir_path, f) for f in readdir(dir_path)])
    
    if isempty(files)
        println(io, "Warning: No txt files found in $dir_path")
        return nothing
    end
    
    # Auto-detect predictor type from first file
    first_file = basename(files[1])
    predictor_prefix = if occursin("rosetta_correlations", first_file)
        "rosetta"
    elseif occursin("foldx_correlations", first_file)
        "foldx"
    elseif occursin("pythia_correlations", first_file)
        "pythia"
    elseif occursin("stabilityoracle_correlations", first_file)
        "stabilityoracle"
    else
        println(io, "Warning: Cannot detect predictor type from filename: $first_file")
        return nothing
    end
    
    println(io, "Detected predictor: $predictor_prefix")
    
    # Define normalization types based on detected predictor
    normalization_types = [
        ("$(predictor_prefix)_correlations", "No normalization"),
        ("n_$(predictor_prefix)_correlations", "Normalization (n)"),
        ("n_n2_$(predictor_prefix)_correlations", "Normalization (n_n2)")
    ]
    
    results = Dict{String, Any}()
    
    println(io, "="^80)
    println(io)
    
    for (prefix, type_name) in normalization_types
        println(io, "Type: $type_name")
        println(io, "-"^80)
        
        # Filter files for this type and target rho
        type_files = filter(f -> startswith(basename(f), prefix) && 
                                 (occursin("rho_$(target_rho).0.txt", basename(f)) ||
                                  occursin("rho_$(target_rho).txt", basename(f))), 
                            files)
        
        if isempty(type_files)
            println(io, "  No files found for rho=$target_rho")
            println(io)
            continue
        end
        
        # Should only be one file for fixed rho
        file = type_files[1]
        
        # Read file content
        content = read(file, String)
        lines = split(content, '\n')
        
        # Extract base correlations
        base_correlations = _extract_base_correlations(lines)
        
        # Extract all combined correlations (for all A values, excluding A=0)
        combined_results = []
        for line in lines
            # Match combined correlation lines
            m = match(r"A=([-\d.]+):\s*Self=([-\d.]+).*Nearby=([-\d.]+)", line)
            if m !== nothing
                A_val = parse(Float64, m.captures[1])
                # Skip A=0 as it's just the original predictor result
                if A_val == 0.0
                    continue
                end
                self_pcc = parse(Float64, m.captures[2])
                nearby_pcc = parse(Float64, m.captures[3])
                push!(combined_results, (A=A_val, Self=self_pcc, Nearby=nearby_pcc))
            end
        end
        
        # Display base correlations
        if !isempty(base_correlations)
            println(io, "\nBase Correlations:")
            # Predictor vs Experimental
            for key in sort(collect(keys(base_correlations)))
                if occursin("_vs_exp", key) && !occursin("ddS", key)
                    predictor_raw = replace(key, "_vs_exp" => "")
                    # Capitalize predictor name properly for display
                    predictor_display = if lowercase(predictor_raw) == "foldx"
                        "FoldX"
                    elseif lowercase(predictor_raw) == "stabilityoracle"
                        "StabilityOracle"
                    else
                        uppercasefirst(lowercase(predictor_raw))
                    end
                    @printf(io, "  %s vs Experimental ddG: %.4f\n", 
                            predictor_display, base_correlations[key])
                end
            end
            # ddS correlations
            if haskey(base_correlations, "self_ddS_vs_exp")
                @printf(io, "  Self ddS vs Experimental ddG: %.4f\n", 
                        base_correlations["self_ddS_vs_exp"])
            end
            if haskey(base_correlations, "nearby_ddS_vs_exp")
                @printf(io, "  Nearby ddS vs Experimental ddG: %.4f\n", 
                        base_correlations["nearby_ddS_vs_exp"])
            end
            if haskey(base_correlations, "self_ddS_vs_diff")
                @printf(io, "  Self ddS vs ddG_diff: %.4f\n", 
                        base_correlations["self_ddS_vs_diff"])
            end
            if haskey(base_correlations, "nearby_ddS_vs_diff")
                @printf(io, "  Nearby ddS vs ddG_diff: %.4f\n", 
                        base_correlations["nearby_ddS_vs_diff"])
            end
        end
        
        # Display combined results - only show best
        if !isempty(combined_results)
            # Find best results
            best_self = maximum(r.Self for r in combined_results)
            best_self_A = [r.A for r in combined_results if r.Self == best_self][1]
            best_nearby = maximum(r.Nearby for r in combined_results)
            best_nearby_A = [r.A for r in combined_results if r.Nearby == best_nearby][1]
            
            println(io, "\nBest Combined Results (Predictor + A*ddS):")
            @printf(io, "  Best Self:   PCC=%.4f, rho=%.1f, A=%.3f\n", best_self, target_rho, best_self_A)
            @printf(io, "  Best Nearby: PCC=%.4f, rho=%.1f, A=%.3f\n", best_nearby, target_rho, best_nearby_A)
            
            results[type_name] = Dict(
                "best_self_pcc" => best_self,
                "best_self_A" => best_self_A,
                "best_nearby_pcc" => best_nearby,
                "best_nearby_A" => best_nearby_A,
                "rho" => target_rho
            )
        else
            println(io, "\nBest Combined Results (Predictor + A*ddS):")
            println(io, "  No combined results found")
        end
        
        println(io, "\n" * "="^80)
        println(io)
    end
    
    return results
end

"""
    load_pcc_results(results_file::String; default_rho=14.0)

Load PCC analysis results from a file.

# Arguments
- `results_file`: Path to the PCC analysis results file
- `default_rho`: Default rho value to use if not found in results

# Returns
- Dict mapping normalization type to named tuple with best parameters
"""
function load_pcc_results(results_file::String)
    if !isfile(results_file)
        error("PCC analysis results file not found: $results_file. Run correlation analysis first.")
    end
    
    println("Reading PCC results from: $results_file")
    
    results = Dict{String, NamedTuple}()
    current_norm = nothing
    
    for line in eachline(results_file)
        if occursin("Type:", line)
            if occursin("No normalization", line) 
                current_norm = "No normalization"
            elseif occursin("Normalization (n_n2)", line) 
                current_norm = "Normalization (n_n2)"
            elseif occursin("Normalization (n)", line) 
                current_norm = "Normalization (n)"
            end
        elseif current_norm !== nothing
            if occursin("Best Self:", line)
                m = match(r"PCC=([-\d.]+).*rho=([-\d.]+).*A=([-\d.]+)", line)
                if m !== nothing
                    curr = get(results, current_norm, (self=0.0, nearby=0.0, pcc_self=0.0, pcc_nearby=0.0, rho_self=0.0, rho_nearby=0.0))
                    results[current_norm] = merge(curr, (self=parse(Float64, m[3]), pcc_self=parse(Float64, m[1]), rho_self=parse(Float64, m[2])))
                end
            elseif occursin("Best Nearby:", line)
                m = match(r"PCC=([-\d.]+).*rho=([-\d.]+).*A=([-\d.]+)", line)
                if m !== nothing
                    curr = get(results, current_norm, (self=0.0, nearby=0.0, pcc_self=0.0, pcc_nearby=0.0, rho_self=0.0, rho_nearby=0.0))
                    results[current_norm] = merge(curr, (nearby=parse(Float64, m[3]), pcc_nearby=parse(Float64, m[1]), rho_nearby=parse(Float64, m[2])))
                end
            end
        end
    end
    
    return results
end

"""
    read_correlation_data(filepath::String, correlation_type::String)

Read correlation data from a correlation results file.

# Arguments
- `filepath`: Path to correlation results file
- `correlation_type`: "self" or "nearby"

# Returns
- Tuple of (A_values, pcc_values) sorted by A values
"""
function read_correlation_data(filepath::String, correlation_type::String)
    A_values = Float64[]
    pcc_values = Float64[]
    
    try
        lines = readlines(filepath)
        
        # Find the section with combined correlations
        in_combined_section = false
        
        for line in lines
            if startswith(line, "Combined correlations")
                in_combined_section = true
                continue
            end
            
            if in_combined_section && startswith(line, "A=")
                first_colon_idx = findfirst(':', line)
                if first_colon_idx !== nothing
                    # Split line into A-value part and correlation part based on the first colon
                    a_value_part = line[1:first_colon_idx-1]
                    corr_part = line[first_colon_idx+1:end]

                    # Extract A value
                    A_str = strip(split(a_value_part, "=")[2])
                    A_val = parse(Float64, A_str)
                    
                    # Now, extract correlation values from the rest of the string (corr_part)
                    if correlation_type == "self"
                        if occursin("Self=", corr_part)
                            self_match = match(r"Self=([+-]?\d*\.?\d+)", corr_part)
                            if self_match !== nothing
                                pcc_val = parse(Float64, self_match.captures[1])
                                push!(A_values, A_val)
                                push!(pcc_values, pcc_val)
                            end
                        end
                    elseif correlation_type == "nearby"
                        if occursin("Nearby=", corr_part)
                            nearby_match = match(r"Nearby=([+-]?\d*\.?\d+)", corr_part)
                            if nearby_match !== nothing
                                pcc_val = parse(Float64, nearby_match.captures[1])
                                push!(A_values, A_val)
                                push!(pcc_values, pcc_val)
                            end
                        end
                    end
                end
            end
        end
    catch e
        println("Error reading correlation file $filepath: $e")
    end
    
    # Sort by A values
    if !isempty(A_values)
        sorted_indices = sortperm(A_values)
        A_values = A_values[sorted_indices]
        pcc_values = pcc_values[sorted_indices]
    end
    
    return A_values, pcc_values
end

"""
    find_best_pcc_parameters(protein_name::String, base_dir::String; threshold::String="70")

Find the best PCC parameters across all predictors for a protein.

# Arguments
- `protein_name`: Name of the protein system
- `base_dir`: Base directory containing the analysis results
- `threshold`: Threshold value (default: "70")

# Returns
- Dict containing best parameters for each predictor
"""
function find_best_pcc_parameters(protein_name::String, base_dir::String; threshold::String="70")
    protein_dir = joinpath(base_dir, protein_name, "figs")
    methods = ["entropy", "entropy_foldx", "entropy_pythia", "entropy_stabilityoracle"]
    
    all_results = Dict{String, Any}()
    
    for method in methods
        dir_path = joinpath(protein_dir, method, threshold)
        if isdir(dir_path)
            results = analyze_pcc_files(dir_path, io=devnull)
            if results !== nothing
                all_results[method] = results
            end
        end
    end
    
    return all_results
end

"""
    save_pcc_results(results::Dict, output_file::String, protein_name::String)

Save PCC analysis results to a file.
"""
function save_pcc_results(results::Dict, output_file::String, protein_name::String)
    mkpath(dirname(output_file))
    
    open(output_file, "w") do io
        println(io, "="^80)
        println(io, "PCC Analysis Results")
        println(io, "Protein: $protein_name")
        println(io, "Date: $(Dates.format(Dates.now(), "yyyy-mm-dd HH:MM:SS"))")
        println(io, "="^80)
        println(io)
        
        for (method, method_results) in results
            println(io, "\n" * "█"^80)
            println(io, "METHOD: $method")
            println(io, "█"^80)
            println(io)
            
            for (norm_type, vals) in method_results
                println(io, "Type: $norm_type")
                println(io, "-"^80)
                
                if haskey(vals, "best_self_pcc") && vals["best_self_pcc"] > -Inf
                    @printf(io, "Best Self:   PCC=%.4f, rho=%.1f, A=%.3f\n", 
                            vals["best_self_pcc"], vals["best_self_rho"], vals["best_self_A"])
                end
                if haskey(vals, "best_nearby_pcc") && vals["best_nearby_pcc"] > -Inf
                    @printf(io, "Best Nearby: PCC=%.4f, rho=%.1f, A=%.3f\n", 
                            vals["best_nearby_pcc"], vals["best_nearby_rho"], vals["best_nearby_A"])
                end
                println(io)
            end
        end
    end
    
    println("Results saved to: $output_file")
end

# Helper function to extract base correlations from lines
function _extract_base_correlations(lines)
    base_correlations = Dict{String, Float64}()
    
    for line in lines
        # Match predictor vs Experimental ddG (case-insensitive)
        m = match(r"(Rosetta|Foldx|Pythia|Stabilityoracle) vs Experimental ddG:\s*([-\d.]+)"i, line)
        if m !== nothing
            predictor = m.captures[1]
            pcc_val = parse(Float64, m.captures[2])
            base_correlations["$(predictor)_vs_exp"] = pcc_val
        end
        
        # Match Self ddS vs Experimental ddG
        m = match(r"Self ddS vs Experimental ddG:\s*([-\d.]+)", line)
        if m !== nothing
            pcc_val = parse(Float64, m.captures[1])
            base_correlations["self_ddS_vs_exp"] = pcc_val
        end
        
        # Match Nearby ddS vs Experimental ddG
        m = match(r"Nearby ddS vs Experimental ddG:\s*([-\d.]+)", line)
        if m !== nothing
            pcc_val = parse(Float64, m.captures[1])
            base_correlations["nearby_ddS_vs_exp"] = pcc_val
        end
        
        # Match Self ddS vs ddG_diff
        m = match(r"Self ddS vs ddG_diff:\s*([-\d.]+)", line)
        if m !== nothing
            pcc_val = parse(Float64, m.captures[1])
            base_correlations["self_ddS_vs_diff"] = pcc_val
        end
        
        # Match Nearby ddS vs ddG_diff
        m = match(r"Nearby ddS vs ddG_diff:\s*([-\d.]+)", line)
        if m !== nothing
            pcc_val = parse(Float64, m.captures[1])
            base_correlations["nearby_ddS_vs_diff"] = pcc_val
        end
    end
    
    return base_correlations
end
