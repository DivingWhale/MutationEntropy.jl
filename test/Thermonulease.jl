using MutationEntropy: phi, gamma
using MutationEntropy
using Test

@testset "msf" begin
    eta1, gamma1 = 1, 1
    eta2, gamma2 = 2, 2
    coordinates = read_coordinates(pkgdir(MutationEntropy, "data", "5XJH.pdb"))
    Γ1 = gamma(coordinates, eta1, gamma1)
    Γ2 = gamma(coordinates, eta2, gamma2)
    @test all(>(0), computed_msf)
end

# Γ = Γ1 + Γ2
# findall(x -> abs(x) < 1e-6, Γ[10, :])

pae = MutationEntropy.read_pae("l190a")