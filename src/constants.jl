#-----------------------------------------------------------------------------
# Amino Acid Size Parameter Mapping (σ_i)
#-----------------------------------------------------------------------------
# These values represent the dimensionless "effective radius" for each
# amino acid residue (Backbone + Side Chain). They are derived by taking the
# cube root of the Total Size Score, which is calculated from the sum of
# all constituent Martini 3 beads. This approach ensures that Glycine has a
# non-zero base size, making the model physically robust.
#
# Formula: σ_i = (Total Score)^(1/3)
#-----------------------------------------------------------------------------

const RESIDUE_SIGMA_MAP = Dict{Char, Float64}(
    # Scaled values (average is 1.0)
    'G' => 0.720,  # Smallest
    'A' => 0.854,
    'S' => 0.907,
    'C' => 0.907,
    'V' => 0.907,
    'P' => 0.907,
    'T' => 0.955,
    'L' => 0.955,
    'I' => 0.955,
    'D' => 0.955,
    'N' => 0.955,
    'M' => 0.998,
    'E' => 0.998,
    'Q' => 0.998,
    'H' => 1.075,
    'K' => 1.075,
    'R' => 1.075,
    'F' => 1.110,
    'Y' => 1.143,
    'W' => 1.203  # Largest
)