function extract_parts(genetic_string::String)
    m= match(r"(\d+)([A-Za-z]*)", genetic_string)
    return m !== nothing ? (parse(Int, m.captures[1]), m.captures[2]) : (nothing, nothing)
end

parse_mutation_position(mutation::AbstractString)::Int = parse(Int, match(r"\d+", mutation).match)

"""
    parse_mutation_details(mutation_str::String)

Parse mutation string (e.g., "a176g") to extract:
- from_aa: wild-type amino acid (single letter, uppercase)
- position: residue position number
- to_aa: mutant amino acid (single letter, uppercase)

Returns (from_aa, position, to_aa)

# Examples
```julia
julia> parse_mutation_details("a176g")
("A", 176, "G")

julia> parse_mutation_details("V42L")
("V", 42, "L")
```
"""
function parse_mutation_details(mutation_str::String)
    m = match(r"^([a-zA-Z]+)(\d+)([a-zA-Z]+)$", mutation_str)
    if m === nothing
        error("Invalid mutation string format: $mutation_str. Expected format like 'a176g' or 'V42L'")
    end
    
    from_aa = uppercase(m.captures[1])
    position = parse(Int, m.captures[2])
    to_aa = uppercase(m.captures[3])
    
    return (from_aa, position, to_aa)
end

"""
    single_to_three_letter(aa::String)

Convert single-letter amino acid code to three-letter code.

# Examples
```julia
julia> single_to_three_letter("A")
"ALA"

julia> single_to_three_letter("w")
"TRP"
```
"""
function single_to_three_letter(aa::String)
    aa_map = Dict(
        "A" => "ALA", "C" => "CYS", "D" => "ASP", "E" => "GLU",
        "F" => "PHE", "G" => "GLY", "H" => "HIS", "I" => "ILE",
        "K" => "LYS", "L" => "LEU", "M" => "MET", "N" => "ASN",
        "P" => "PRO", "Q" => "GLN", "R" => "ARG", "S" => "SER",
        "T" => "THR", "V" => "VAL", "W" => "TRP", "Y" => "TYR"
    )
    return get(aa_map, uppercase(aa), aa)  # Return original if not found
end

# Amino acid size-based ρ calculation
# Based on relative amino acid sizes from Grantham (1974) and other biochemical references
"""
    AMINO_ACID_RHO

Dictionary mapping three-letter amino acid codes to size-based ρ values.
Smaller amino acids have larger ρ values, larger amino acids have smaller ρ values.
Based on relative amino acid sizes from Grantham (1974) and biochemical properties.
"""
const AMINO_ACID_RHO = Dict(
    "GLY" => 12.00,  # Size: 1.00 (smallest)
    "ALA" => 11.45,  # Size: 1.67
    "SER" => 11.18,  # Size: 2.00
    "CYS" => 11.18,  # Size: 2.00
    "VAL" => 11.18,  # Size: 2.00
    "PRO" => 11.18,  # Size: 2.00
    "THR" => 10.91,  # Size: 2.33
    "LEU" => 10.91,  # Size: 2.33
    "ILE" => 10.91,  # Size: 2.33
    "ASP" => 10.91,  # Size: 2.33
    "ASN" => 10.91,  # Size: 2.33
    "MET" => 10.64,  # Size: 2.67
    "GLU" => 10.64,  # Size: 2.67
    "GLN" => 10.64,  # Size: 2.67
    "HIS" => 10.09,  # Size: 3.33
    "LYS" => 10.09,  # Size: 3.33
    "ARG" => 10.09,  # Size: 3.33
    "PHE" => 9.82,   # Size: 3.67
    "TYR" => 9.55,   # Size: 4.00
    "TRP" => 9.00    # Size: 4.67 (largest)
)

"""
    get_residue_rho(residue_name::String)

Get the size-based ρ value for a given amino acid residue.
Accepts both single-letter and three-letter amino acid codes.

# Arguments
- `residue_name::String`: Amino acid code (single or three letter)

# Returns
- `Float64`: The ρ value for the amino acid (default: 10.5 for unknown)

# Examples
```julia
julia> get_residue_rho("ALA")
11.45

julia> get_residue_rho("A")
11.45

julia> get_residue_rho("TRP")
9.0
```
"""
function get_residue_rho(residue_name::AbstractString)
    aa_code = uppercase(String(residue_name))
    
    # If single letter, convert to three letter
    if length(aa_code) == 1
        aa_code = single_to_three_letter(aa_code)
    end
    
    return get(AMINO_ACID_RHO, aa_code, 10.5)  # Default to mid-range value
end

function read_mutations_from_file(task_file_path::String)::Vector{String}
    return [strip(line) for line in eachline(task_file_path) if !isempty(strip(line))]
end

function get_experimental_ddg(ddG_exp::DataFrame, position::Int, mutation_residue::String)::Union{Float64, Nothing}
    matching_rows = (ddG_exp.position .== position) .& (ddG_exp.mutation .== mutation_residue)
    
    if !any(matching_rows)
        return nothing
    end
    
    return mean(ddG_exp[matching_rows, :ddG])
end

function _ensure_pdb_from_cif_if_needed(subfolder::String, model_pdb_path::String)
    if isfile(model_pdb_path)
        return true # PDB already exists
    end

    # Check if CIF file exists as a source for conversion
    cif_file = joinpath(subfolder, "model.cif")
    if !isfile(cif_file)
        
        error("Neither PDB nor CIF file found in $subfolder") 
    end

    # Convert CIF to PDB using Python script with conda environment
    script_path = abspath(joinpath(@__DIR__, "convert.py"))
    subfolder_abs = abspath(subfolder)
    conda_env = "bio" # As used in the original code
    
    # Find conda executable
    conda_paths = [
        "/opt/miniconda3/bin/conda",
        "/opt/homebrew/bin/conda",
        "/usr/local/bin/conda", 
        "/Users/sam/miniconda3/bin/conda",
        expanduser("~/miniconda3/bin/conda"),
        expanduser("~/anaconda3/bin/conda")
    ]
    
    conda_exe = nothing
    for path in conda_paths
        if isfile(path)
            conda_exe = path
            break
        end
    end
    
    if conda_exe === nothing
        # Fallback: try to find conda in PATH
        try
            conda_exe = strip(read(`which conda`, String))
        catch
            error("Could not find conda executable. Please ensure conda is installed and accessible.")
        end
    end

    cmd = `$conda_exe run -n $conda_env python $script_path $subfolder_abs`

    try
        run(cmd)
    catch e
        println("Error during conversion script execution: $e") 
        error("Failed to convert CIF to PDB. Check conda environment '$conda_env', Python script '$script_path', and related files.")
    end
    
    if !isfile(model_pdb_path)
        error("Conversion attempted but PDB file not created: $model_pdb_path")
    end
    return true
end

"""
    get_low_plddt_residues(mutation::String, round::Int, datadir::String; threshold::Float64=90.0)

Identify residues with low pLDDT scores.
"""
function get_low_plddt_residues(mutation::String, round_val::Int, datadir::String; threshold::Float64=90.0) # Renamed round
    # Construct the base directory path
    base_dir = joinpath(datadir, "$(mutation)_$(round_val)")
    
    # Use glob with relative pattern from the base directory
    subfolder_pattern = "seed-*_sample-0"
    
    if !isdir(base_dir)
        error("Base directory not found: $base_dir")
    end
    
    subfolders = glob(subfolder_pattern, base_dir)
    # Convert relative paths to absolute paths
    subfolders = [joinpath(base_dir, sf) for sf in subfolders]
    
    if isempty(subfolders)
        error("No matching folder found: $subfolder_pattern in $base_dir")
    end
    
    subfolder = subfolders[1]
    model_pdb = joinpath(subfolder, "model.pdb")
    # temp_pdb = "$(mutation)_temp_pl.pdb" # This variable was unused in the original.
    
    _ensure_pdb_from_cif_if_needed(subfolder, model_pdb) # Handles PDB creation

    structure = read(model_pdb, PDBFormat)
    
    confidence_file = joinpath(subfolder, "confidences.json")
    if !isfile(confidence_file)
        error("Confidence file not found: $confidence_file")
    end
    
    confidence_data = JSON.parsefile(confidence_file)
    # Ensure atom_plddts is Vector{Float64} for `mean` function compatibility
    atom_plddts = convert(Vector{Float64}, confidence_data["atom_plddts"])
    
    residue_plddts = Dict{Int, Float64}()
    atom_idx_tracker = 1 # Use a different name to avoid confusion with atom_index if it's a global/module var
    
    for model_obj in structure # model is a type in BioStructures
        for chain_obj in model_obj # chain is a type in BioStructures
            for residue_obj in chain_obj # residue is a type in BioStructures
                # Original parsing logic, will error if resid is not parsable to Int
                res_id = parse(Int, BioStructures.resid(residue_obj))
                
                # Using collectatoms for consistency, original had collect(atoms(residue))
                current_residue_atoms = collectatoms(residue_obj)
                atom_count = length(current_residue_atoms)
                
                if atom_count > 0 && !haskey(residue_plddts, res_id)
                    # Ensure atom_idx_tracker and atom_count do not exceed atom_plddts bounds
                    if atom_idx_tracker + atom_count - 1 <= lastindex(atom_plddts)
                        plddt_scores_for_residue = view(atom_plddts, atom_idx_tracker:(atom_idx_tracker + atom_count - 1))
                        residue_plddts[res_id] = mean(plddt_scores_for_residue)
                    else
                        @warn "Atom index out of bounds for residue $res_id (mutation $mutation, round $round_val). Atom pLDDTs length: $(length(atom_plddts)), trying to access up to $(atom_idx_tracker + atom_count - 1). Skipping pLDDT for this residue."
                    end
                elseif atom_count == 0 && !haskey(residue_plddts, res_id) # Handle residues with no atoms if they occur
                     @warn "Residue $res_id (mutation $mutation, round $round_val) has no atoms. Skipping pLDDT calculation."
                end
                
                # Crucially, advance atom_idx_tracker by atom_count for every residue processed,
                # regardless of whether its pLDDT was stored, to maintain sync with atom_plddts list.
                # This fixes a potential bug in the original atom_index advancement.
                atom_idx_tracker += atom_count
            end
        end
    end
    
    low_plddt_res_ids = [res_id for (res_id, plddt) in residue_plddts if plddt < threshold]
    
    return sort(low_plddt_res_ids)
end

"""
    find_residues_within_distance(residue_id::Int, d_matrix::Matrix{Float64}; distance::Float64=13.0)

Find all residues within a specific distance of a target residue using matrix indices.

Args:
- residue_id: Matrix index of target residue (1-based, for truncated matrix)
- d_matrix: Distance matrix (truncated)
- distance: Distance threshold in Angstroms

Returns:
- Vector of matrix indices for residues within distance threshold

Note: Returns indices including self-distance (0.0) and handles matrix bounds checking.
Subsequent processing typically filters out the target residue itself.
"""
function find_residues_within_distance(residue_id::Int, d_matrix::Matrix{Float64}; distance::Float64=13.0)
    # if !(1 <= residue_id <= size(d_matrix, 1))
    #     @warn "Residue ID $residue_id is out of bounds for distance matrix row access (max rows: $(size(d_matrix,1))). Returning empty list."
    #     return Int[]
    # end
    # Original logic: findall(x -> x < distance, d_matrix[residue_id, :])
    # This includes self-distance (0.0) if distance > 0.
    # Subsequent code (calculate_strain, collect_strains) handles `if residue == site; continue; end`.
    nearby_indices = findall(x -> x < distance, d_matrix[residue_id, :])
    
    return sort(nearby_indices) # Returns sorted Vector{Int}
end

"""
    read_coordinates(filename)

Reads CA atom coordinates from a PDB file.
"""
function read_coordinates(filename::String)::Vector{Vector{Float64}}
    structure = read(filename, PDBFormat)
    ca_atoms = collectatoms(structure, sel"name CA")
    return getfield.(ca_atoms, :coords)
end

"""
    read_xvg(filename)

Parses a generic .xvg file, returning a vector of (Int, Float64) tuples.
"""
function read_xvg(filename::String)::Vector{Tuple{Int, Float64}}
    results = Tuple{Int, Float64}[]
    for line in eachline(filename)
        if startswith(line, ('#', '@')) || isempty(line)
            continue
        end
        parts = split(line)
        if length(parts) >= 2
            push!(results, (parse(Int, parts[1]), parse(Float64, parts[2])))
        end
    end
    return results
end