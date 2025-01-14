module MutationEntropy

using BioStructures, LinearAlgebra
using CairoMakie

export read_coordinates, read_xvg, plot_msf

include("entropy.jl")
include("visualization.jl")

end
