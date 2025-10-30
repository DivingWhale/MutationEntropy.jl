"""
    MutationData

Contains all matrices and identifiers for a single mutation analysis.
Includes sequence information for sigma-weighted entropy calculations.
"""
struct MutationData
    wt_pae::Vector{Matrix{Float64}}
    mutant_pae::Vector{Matrix{Float64}}
    wt_dist::Vector{Matrix{Float64}}
    mutant_dist::Vector{Matrix{Float64}}
    mutation::String
    wt_sequence::String           # WT sequence (one-letter codes)
    mutation_residue::Char        # Mutated residue (one-letter code)
end

# Constructor without sequence information for backward compatibility
function MutationData(
    wt_pae::Vector{Matrix{Float64}},
    mutant_pae::Vector{Matrix{Float64}},
    wt_dist::Vector{Matrix{Float64}},
    mutant_dist::Vector{Matrix{Float64}},
    mutation::String
)
    return MutationData(wt_pae, mutant_pae, wt_dist, mutant_dist, mutation, "", ' ')
end

"""
    EntropyParams

Parameters for entropy calculation with simplified indexing system.

Fields:
- `position::Int`: Biological residue number for calculation (e.g., 88, 89, 90...)
- `rho::Float64`: PAE weighting parameter  
- `α::Float64`: Distance weighting parameter
- `matrix_start::Int`: First biological residue number in the truncated matrix (e.g., 88)
                       Used to convert: matrix_idx = position - matrix_start + 1
- `filter_low_plddt::Bool`: Whether to filter low pLDDT residues
- `plddt_threshold::Float64`: pLDDT threshold for filtering
- `data_dir::String`: Directory containing pLDDT data

Note: For backward compatibility with legacy code, matrix_start can be provided
      as (offset + 1) where offset was the old parameter.
"""
struct EntropyParams
    position::Int
    rho::Float64
    α::Float64
    matrix_start::Int  # Renamed from offset for clarity
    filter_low_plddt::Bool
    plddt_threshold::Float64
    data_dir::String
end
