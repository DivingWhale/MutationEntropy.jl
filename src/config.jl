"""
Protein System Configuration Module

This module provides a standardized way to configure protein analysis systems.
Each protein system can define its configuration in a simple TOML file or Julia struct.
"""

using TOML

export ProteinConfig, PredictorConfig, load_protein_config, get_predictor_config
export DEFAULT_ANALYSIS_PARAMS

"""
Configuration for a single predictor (e.g., Rosetta, FoldX, Pythia)
"""
struct PredictorConfig
    name::String
    read_func::String
    path::String
    output_dir::String
    pdb_id::Union{String, Nothing}
    offset::Union{Int, Nothing}
    entropy_offset::Union{Int, Nothing}  # Predictor-specific offset to convert entropy numbering to predictor numbering
end

"""
Configuration for an entire protein analysis system
"""
struct ProteinConfig
    name::String
    base_data_dir::String
    project_dir::String
    pae_path::String
    experimental_data_path::String
    figures_dir::String
    reference_pdb::Union{String, Nothing}
    pdb_offset::Int
    entropy_offset::Int  # Offset to convert entropy mutant numbering to experimental/PDB numbering
    negate_exp_ddg::Bool  # If true, negate experimental ddG values (for ThermoMutDB sign convention)
    predictors::Dict{String, PredictorConfig}
    # Analysis parameters
    target_alphas::Vector{Float64}
    target_rhos::Vector{Float64}
    threshold::String
    a_values::Vector{Float64}
    normalization_configs::Vector{Tuple{Bool, Bool}}
end

# Default analysis parameters
const DEFAULT_ANALYSIS_PARAMS = (
    target_alphas = [0.0],
    target_rhos = [8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0],
    threshold = "70",
    a_values = vcat(-1:0.1:-0.2, -0.1:0.01:-0.02, -0.01:0.001:0),
    normalization_configs = [(false, false), (true, false), (true, true)]
)

"""
    load_protein_config(config_path::String) -> ProteinConfig

Load protein configuration from a TOML file.

Example TOML structure:

    [protein]
    name = "1PGA"
    base_data_dir = "/path/to/data"
    project_dir = "/path/to/project"
    reference_pdb = "/path/to/ref.pdb"
    pdb_offset = -226

    [predictors.rosetta]
    read_func = "read_ddg_rosetta"
    path = "\${base_data_dir}/Rosetta"
    output_dir = "entropy"

    [predictors.foldx]
    read_func = "read_ddg_foldx"
    path = "\${base_data_dir}/FoldX"
    output_dir = "entropy_foldx"
    pdb_id = "1PGA"
"""
function load_protein_config(config_path::String)::ProteinConfig
    if !isfile(config_path)
        error("Configuration file not found: $config_path")
    end
    
    config = TOML.parsefile(config_path)
    
    # Parse protein section
    protein = config["protein"]
    base_data_dir = protein["base_data_dir"]
    project_dir = get(protein, "project_dir", dirname(config_path))
    
    # Parse predictors
    predictors = Dict{String, PredictorConfig}()
    if haskey(config, "predictors")
        for (name, pred_config) in config["predictors"]
            # Replace ${base_data_dir} placeholder
            path = replace(get(pred_config, "path", ""), "\${base_data_dir}" => base_data_dir)
            
            predictors[name] = PredictorConfig(
                name,
                pred_config["read_func"],
                path,
                get(pred_config, "output_dir", "entropy_$name"),
                get(pred_config, "pdb_id", nothing),
                get(pred_config, "offset", nothing),
                get(pred_config, "entropy_offset", nothing)  # Predictor-specific entropy offset
            )
        end
    end
    
    # Parse analysis parameters with defaults
    analysis = get(config, "analysis", Dict())
    target_alphas = get(analysis, "target_alphas", DEFAULT_ANALYSIS_PARAMS.target_alphas)
    target_rhos = get(analysis, "target_rhos", DEFAULT_ANALYSIS_PARAMS.target_rhos)
    threshold = get(analysis, "threshold", DEFAULT_ANALYSIS_PARAMS.threshold)
    
    # Parse a_values - supports string range like "-10:0.1:0" or array
    a_values_raw = get(analysis, "a_values", nothing)
    a_values = if a_values_raw === nothing
        DEFAULT_ANALYSIS_PARAMS.a_values
    elseif a_values_raw isa String
        # Parse range string like "-10:0.1:0" or "vcat(-10:0.1:0)"
        range_str = replace(a_values_raw, r"^vcat\((.+)\)$" => s"\1")
        parts = split(range_str, ":")
        if length(parts) == 3
            collect(range(parse(Float64, parts[1]), parse(Float64, parts[3]), step=parse(Float64, parts[2])))
        elseif length(parts) == 2
            collect(range(parse(Float64, parts[1]), parse(Float64, parts[2])))
        else
            error("Invalid a_values range format: $a_values_raw")
        end
    else
        Float64.(a_values_raw)
    end
    
    return ProteinConfig(
        protein["name"],
        base_data_dir,
        project_dir,
        get(protein, "pae_path", joinpath(base_data_dir, "af3")),
        get(protein, "experimental_data_path", joinpath(base_data_dir, "$(protein["name"]).csv")),
        get(protein, "figures_dir", joinpath(project_dir, "figs")),
        get(protein, "reference_pdb", nothing),
        get(protein, "pdb_offset", 0),
        get(protein, "entropy_offset", 0),  # Offset for entropy mutant numbering
        get(protein, "negate_exp_ddg", false),  # Negate experimental ddG (for ThermoMutDB)
        predictors,
        Float64.(target_alphas),
        Float64.(target_rhos),
        string(threshold),
        a_values,
        DEFAULT_ANALYSIS_PARAMS.normalization_configs
    )
end

"""
    create_protein_config(; kwargs...) -> ProteinConfig

Create a ProteinConfig programmatically from keyword arguments.
"""
function create_protein_config(;
    name::String,
    base_data_dir::String,
    project_dir::String = "",
    pae_path::String = "",
    experimental_data_path::String = "",
    figures_dir::String = "",
    reference_pdb::Union{String, Nothing} = nothing,
    pdb_offset::Int = 0,
    entropy_offset::Int = 0,
    negate_exp_ddg::Bool = false,
    predictors::Dict = Dict(),
    target_alphas::Vector{Float64} = DEFAULT_ANALYSIS_PARAMS.target_alphas,
    target_rhos::Vector{Float64} = DEFAULT_ANALYSIS_PARAMS.target_rhos,
    threshold::String = DEFAULT_ANALYSIS_PARAMS.threshold,
    a_values::Vector{Float64} = DEFAULT_ANALYSIS_PARAMS.a_values,
    normalization_configs = DEFAULT_ANALYSIS_PARAMS.normalization_configs
)::ProteinConfig
    
    # Set defaults based on name and base_data_dir
    project_dir = isempty(project_dir) ? joinpath(dirname(base_data_dir), name) : project_dir
    pae_path = isempty(pae_path) ? joinpath(base_data_dir, "af3") : pae_path
    experimental_data_path = isempty(experimental_data_path) ? joinpath(base_data_dir, "$name.csv") : experimental_data_path
    figures_dir = isempty(figures_dir) ? joinpath(project_dir, "figs") : figures_dir
    
    # Convert predictor dict if needed
    pred_configs = Dict{String, PredictorConfig}()
    for (pname, pconfig) in predictors
        if pconfig isa PredictorConfig
            pred_configs[pname] = pconfig
        elseif pconfig isa NamedTuple || pconfig isa Dict
            pred_configs[pname] = PredictorConfig(
                pname,
                pconfig[:read_func],
                pconfig[:path],
                get(pconfig, :output_dir, "entropy_$pname"),
                get(pconfig, :pdb_id, nothing),
                get(pconfig, :offset, nothing),
                get(pconfig, :entropy_offset, nothing)
            )
        end
    end
    
    return ProteinConfig(
        name, base_data_dir, project_dir, pae_path, experimental_data_path,
        figures_dir, reference_pdb, pdb_offset, entropy_offset, negate_exp_ddg, pred_configs,
        target_alphas, target_rhos, threshold, a_values, collect(normalization_configs)
    )
end

"""
    get_predictor_config(config::ProteinConfig, predictor_name::String) -> PredictorConfig

Get the configuration for a specific predictor from a ProteinConfig.
"""
function get_predictor_config(config::ProteinConfig, predictor_name::String)
    if !haskey(config.predictors, predictor_name)
        error("Predictor '$predictor_name' not found. Available: $(keys(config.predictors))")
    end
    return config.predictors[predictor_name]
end

"""
    config_from_legacy(config_local_path::String) -> ProteinConfig

Convert a legacy config_local.jl file to ProteinConfig.
This function includes the legacy config file and extracts the configuration.
"""
function config_from_legacy(config_local_path::String)
    if !isfile(config_local_path)
        error("Config file not found: $config_local_path")
    end
    
    # This creates a temporary module to safely include the config
    config_module = Module()
    Base.include(config_module, config_local_path)
    
    # Extract variables from the included module
    base_data_dir = getfield(config_module, :BASE_DATA_DIR)
    project_dir = getfield(config_module, :PROJECT_DIR)
    
    # Get optional fields with defaults
    pae_path = hasproperty(config_module, :PAE_PATH) ? getfield(config_module, :PAE_PATH) : joinpath(base_data_dir, "af3")
    exp_path = getfield(config_module, :EXPERIMENTAL_DATA_PATH)
    figures_dir = hasproperty(config_module, :FIGURES_DIR) ? getfield(config_module, :FIGURES_DIR) : joinpath(project_dir, "figs")
    
    # Get predictor configs
    legacy_predictors = getfield(config_module, :PREDICTOR_CONFIGS)
    
    predictors = Dict{String, PredictorConfig}()
    for (name, pred) in legacy_predictors
        predictors[name] = PredictorConfig(
            name,
            pred.read_func,
            pred.path,
            pred.output_dir,
            hasproperty(pred, :PDB_ID) ? pred.PDB_ID : nothing,
            hasproperty(pred, :offset) ? pred.offset : nothing,
            hasproperty(pred, :entropy_offset) ? pred.entropy_offset : nothing
        )
    end
    
    # Extract protein name from project dir
    protein_name = basename(project_dir)
    
    return create_protein_config(
        name = protein_name,
        base_data_dir = base_data_dir,
        project_dir = project_dir,
        pae_path = pae_path,
        experimental_data_path = exp_path,
        figures_dir = figures_dir,
        predictors = predictors
    )
end

"""
    generate_config_toml(config::ProteinConfig, output_path::String)

Generate a TOML configuration file from a ProteinConfig struct.
"""
function generate_config_toml(config::ProteinConfig, output_path::String)
    lines = String[]
    
    push!(lines, "# Protein Configuration File")
    push!(lines, "# Generated automatically - modify as needed")
    push!(lines, "")
    push!(lines, "[protein]")
    push!(lines, "name = \"$(config.name)\"")
    push!(lines, "base_data_dir = \"$(config.base_data_dir)\"")
    push!(lines, "project_dir = \"$(config.project_dir)\"")
    push!(lines, "pae_path = \"$(config.pae_path)\"")
    push!(lines, "experimental_data_path = \"$(config.experimental_data_path)\"")
    push!(lines, "figures_dir = \"$(config.figures_dir)\"")
    if config.reference_pdb !== nothing
        push!(lines, "reference_pdb = \"$(config.reference_pdb)\"")
    end
    push!(lines, "pdb_offset = $(config.pdb_offset)")
    push!(lines, "entropy_offset = $(config.entropy_offset)")
    push!(lines, "")
    
    push!(lines, "[analysis]")
    push!(lines, "target_alphas = [$(join(config.target_alphas, ", "))]")
    push!(lines, "target_rhos = [$(join(config.target_rhos, ", "))]")
    push!(lines, "threshold = \"$(config.threshold)\"")
    push!(lines, "")
    
    for (name, pred) in config.predictors
        push!(lines, "[predictors.$name]")
        push!(lines, "read_func = \"$(pred.read_func)\"")
        push!(lines, "path = \"$(pred.path)\"")
        push!(lines, "output_dir = \"$(pred.output_dir)\"")
        if pred.pdb_id !== nothing
            push!(lines, "pdb_id = \"$(pred.pdb_id)\"")
        end
        if pred.offset !== nothing
            push!(lines, "offset = $(pred.offset)")
        end
        if pred.entropy_offset !== nothing
            push!(lines, "entropy_offset = $(pred.entropy_offset)")
        end
        push!(lines, "")
    end
    
    mkpath(dirname(output_path))
    write(output_path, join(lines, "\n"))
    println("Generated config file: $output_path")
end
