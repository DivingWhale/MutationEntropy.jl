# load data
using BioStructures, LinearAlgebra
using Test

phi(r, eta, kappa) = exp(-(r/eta)^kappa)

function gamma(coordinates::AbstractVector{Vector{Float64}}, eta::Int, gamma::Int)
    natoms = length(coordinates)
    Γ = zeros(natoms, natoms)
    for i = 1:natoms, j = i+1:natoms
        Γ[i, j] = -phi(norm(coordinates[i] - coordinates[j]), eta, gamma)
    end
    Γ += Γ'
    for i = 1:natoms
        Γ[i, i] = -sum(view(Γ, :, i))
    end
    return Γ
end

function read_coordinates(filename::String)
    WT_structure = read(filename, PDBFormat)
    CA_atoms = collectatoms(WT_structure, sel"name CA")
    coordinates = getfield.(CA_atoms, :coords)
    return coordinates
end

@testset "gamma" begin
    coordinates = read_coordinates("./5XJH.pdb")
    Γ = gamma(coordinates, 1, 1)
    @test Γ' ≈ Γ
    @test all(x->isapprox(x, 0; atol=1e-8), sum(Γ, dims=1))
end

function read_xvg(filename::String)
    result = Tuple{Int, Float64}[]
    for line in eachline(filename)
        if startswith(line, "#") || startswith(line, "@") || isempty(line)
            continue
        end
        parts = split(line)
        push!(result, (parse(Int, parts[1]), parse(Float64, parts[2])))
    end
    return result
end

function msf(Γ::AbstractMatrix)
    return diag(pinv(Γ))
end

##
eta1, gamma1 = 1, 1
eta2, gamma2 = 2, 2
coordinates = read_coordinates("./5XJH.pdb")
Γ1 = gamma(coordinates, eta1, gamma1)
Γ2 = gamma(coordinates, eta2, gamma2)
computed_msf = msf(Γ1 + Γ2)

##
experimental_data = read_xvg("5xjh_A298_10us_cat12_rmsf_bb.xvg")
rmsf_exp = getindex.(experimental_data, 2)
msf_exp = rmsf_exp.^(-2)

# Plot computed_msf vs msf_exp, with CairoMakie
using CairoMakie
fig = Figure()
ax = Axis(fig[1, 1])
lines!(ax, msf_exp; label="Experimental MSF")
lines!(ax, computed_msf; label="Computed MSF")
# xlabel!(ax, "Experimental MSF")
# ylabel!(ax, "Computed MSF")
fig