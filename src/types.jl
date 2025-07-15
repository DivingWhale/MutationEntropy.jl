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

Parameters for entropy calculation.
"""
struct EntropyParams
    position::Int
    rho::Float64
    α::Float64
    offset::Int
    filter_low_plddt::Bool
    plddt_threshold::Float64
    data_dir::String
end
