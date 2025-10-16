"""
    MutationData

Contains all matrices and identifiers for a single mutation analysis.
"""
struct MutationData
    wt_pae::Vector{Matrix{Float64}}
    mutant_pae::Vector{Matrix{Float64}}
    wt_dist::Vector{Matrix{Float64}}
    mutant_dist::Vector{Matrix{Float64}}
    mutation::String
end

"""
    EntropyParams

Parameters for entropy calculation with dual ρ support for mutant and wild-type.

Fields:
- `position::Int`: Biological residue number for calculation (e.g., 88, 89, 90...)
- `rho_mut::Float64`: PAE weighting parameter for mutant (size-based)
- `rho_wt::Float64`: PAE weighting parameter for wild-type (size-based)
- `α::Float64`: Distance weighting parameter
- `offset::Int`: Offset used to convert position to matrix index (matrix_idx = position - offset)
- `filter_low_plddt::Bool`: Whether to filter low pLDDT residues
- `plddt_threshold::Float64`: pLDDT threshold for filtering
- `data_dir::String`: Directory containing pLDDT data

"""
struct EntropyParams
    position::Int
    rho_mut::Float64
    rho_wt::Float64
    α::Float64
    offset::Int
    filter_low_plddt::Bool
    plddt_threshold::Float64
    data_dir::String
end
