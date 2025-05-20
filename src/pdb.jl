using PDBTools

wt_atoms = read_pdb("../5XJH.pdb")
ca_atoms = select(wt_atoms, "name CA")

dis_matrix = contact_map(ca_atoms; discrete=false, dmax=Inf, gap=3)

dis_matrix.matrix


function replace_missing_with_nan!(matrix::Matrix{Union{Missing, Float32}})::Matrix{Float32}
    """
    Replaces all Missing values in a Matrix{Union{Missing, Float32}} with NaN and returns a new Matrix{Float32}.
    This function modifies the input matrix in place.
  
    Args:
      matrix: The input Matrix{Union{Missing, Float32}}.
  
    Returns:
      A new Matrix{Float32} with Missing values replaced by NaN.
      Note: The function modifies the original matrix in-place.
    """
  
    # Create a new matrix of Float32 to hold the result
    new_matrix = Matrix{Float32}(undef, size(matrix))
  
    for i in eachindex(matrix)
      if ismissing(matrix[i])
        new_matrix[i] = NaN32
      else
        new_matrix[i] = Float32(matrix[i]) # Convert to Float32 explicitly
      end
    end
  
    return new_matrix
  end
