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
    # Total Score = 3
    'G' => 1.442,

    # Total Score = 5
    'A' => 1.710,

    # Total Score = 6
    'S' => 1.817,
    'C' => 1.817,
    'V' => 1.817,
    'P' => 1.817,

    # Total Score = 7
    'T' => 1.913,
    'L' => 1.913,
    'I' => 1.913,
    'D' => 1.913,
    'N' => 1.913,

    # Total Score = 8
    'M' => 2.000,
    'E' => 2.000,
    'Q' => 2.000,

    # Total Score = 10
    'H' => 2.154,
    'K' => 2.154,
    'R' => 2.154,

    # Total Score = 11
    'F' => 2.224,

    # Total Score = 12
    'Y' => 2.289,

    # Total Score = 14
    'W' => 2.410
)
