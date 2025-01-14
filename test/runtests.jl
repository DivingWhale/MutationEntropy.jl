using MutationEntropy: phi, gamma
using MutationEntropy
using Test

@testset "gamma" begin
    coordinates = read_coordinates(pkgdir(MutationEntropy, "data", "5XJH.pdb"))
    Γ = gamma(coordinates, 1, 1)
    @test Γ' ≈ Γ
    @test all(x->isapprox(x, 0; atol=1e-8), sum(Γ, dims=1))
end

@testset "msf" begin
    eta1, gamma1 = 1, 1
    eta2, gamma2 = 2, 2
    coordinates = read_coordinates(pkgdir(MutationEntropy, "data", "5XJH.pdb"))
    Γ1 = gamma(coordinates, eta1, gamma1)
    Γ2 = gamma(coordinates, eta2, gamma2)
    computed_msf = MutationEntropy.msf(Γ1 + Γ2)
    @test all(>(0), computed_msf)
end

@testset "read xvg" begin
    experimental_data = MutationEntropy.read_xvg(pkgdir(MutationEntropy, "data", "5xjh_A298_10us_cat12_rmsf_bb.xvg"))
    rmsf_exp = getindex.(experimental_data, 2)
    msf_exp = rmsf_exp.^(-2)
    @test all(>(0), msf_exp)
end
