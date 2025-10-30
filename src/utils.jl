function extract_parts(genetic_string::String)
    m= match(r"(\d+)([A-Za-z]*)", genetic_string)
    return m !== nothing ? (parse(Int, m.captures[1]), m.captures[2]) : (nothing, nothing)
end

parse_mutation_position(mutation::AbstractString)::Int = parse(Int, match(r"\d+", mutation).match)

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

#-----------------------------------------------------------------------------
# Sequence Extraction and Sigma Calculation Functions
#-----------------------------------------------------------------------------

"""
    extract_sequence_from_structure(structure)

Extract one-letter amino acid sequence from a BioStructures structure.
Uses BioSequences.LongAA for robust sequence extraction.
Returns a string of one-letter amino acid codes.

Note: The sequence starts from the first residue present in the PDB file.
For Thermonuclease, PDB position 6 = biological 88, so sequence[1] corresponds
to biological residue 88.
"""
function extract_sequence_from_structure(structure)::String
    # Get the first model and first chain
    model = first(structure)
    chain = first(model)
    
    # Extract sequence using LongAA with standard selector
    # standardselector filters out non-standard residues
    seq = LongAA(chain, standardselector)
    
    return String(seq)
end

"""
    parse_mutation_residue(mutation::String)

Extract the mutated residue (one-letter code) from a mutation string.
Example: "e75v" -> 'V'
"""
function parse_mutation_residue(mutation::String)::Char
    # Mutation format: original_residue + position + mutated_residue
    # e.g., "e75v" means position 75, E->V
    mutation_lower = lowercase(mutation)
    if length(mutation_lower) < 2
        error("Invalid mutation format: $mutation")
    end
    
    # Last character is the mutated residue
    mutated_residue = uppercase(mutation_lower[end])
    
    return mutated_residue
end

"""
    get_sigma(residue::Char)

Get the sigma value for a given amino acid residue.
Returns NaN for unknown residues.
"""
function get_sigma(residue::Char)::Float64
    residue_upper = uppercase(residue)
    if haskey(RESIDUE_SIGMA_MAP, residue_upper)
        return RESIDUE_SIGMA_MAP[residue_upper]
    else
        @warn "Unknown residue for sigma lookup: $residue. Returning NaN."
        return NaN
    end
end

"""
    calculate_sigma_ij(idx_i::Int, idx_j::Int, sequence::String, mutation_position::Int, mutation_residue::Char, matrix_start::Int)

Calculate the average sigma value for a residue pair.
For mutant calculations, replaces the WT residue at mutation_position with mutation_residue.

Args:
- idx_i, idx_j: Matrix indices (1-based, where 1 = first position in truncated matrix)
- sequence: WT sequence extracted from PDB (already truncated to match matrix)
- mutation_position: Biological residue number of mutation (e.g., 88 for Thermonuclease)
- mutation_residue: Mutated residue (one-letter code), use ' ' for WT
- matrix_start: First biological residue number in matrix (e.g., 88 for Thermonuclease)

Simplified indexing:
- matrix_idx 1 = biological position matrix_start
- sequence[1] = biological position matrix_start
- biological_position = matrix_idx + matrix_start - 1
"""
function calculate_sigma_ij(
    idx_i::Int, 
    idx_j::Int, 
    sequence::String, 
    mutation_position::Int, 
    mutation_residue::Char,
    matrix_start::Int
)::Float64
    # Matrix indices directly correspond to sequence positions
    # because both start from the same truncation point
    seq_pos_i = idx_i
    seq_pos_j = idx_j
    
    # Check bounds
    if seq_pos_i < 1 || seq_pos_i > length(sequence) || seq_pos_j < 1 || seq_pos_j > length(sequence)
        @warn "Sequence position out of bounds: i=$seq_pos_i, j=$seq_pos_j, sequence length=$(length(sequence))"
        return NaN
    end
    
    # Get residues from WT sequence
    residue_i = sequence[seq_pos_i]
    residue_j = sequence[seq_pos_j]
    
    # If this is a mutant calculation (mutation_residue != ' '), 
    # replace the residue at the mutation site
    if mutation_residue != ' '
        # Convert matrix indices to biological positions
        bio_pos_i = idx_i + matrix_start - 1
        bio_pos_j = idx_j + matrix_start - 1
        
        if bio_pos_i == mutation_position
            residue_i = mutation_residue
        end
        if bio_pos_j == mutation_position
            residue_j = mutation_residue
        end
    end
    
    sigma_i = get_sigma(residue_i)
    sigma_j = get_sigma(residue_j)
    
    if isnan(sigma_i) || isnan(sigma_j)
        return NaN
    end
    
    return (sigma_i + sigma_j) / 2.0
end
