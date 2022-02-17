@testset "NVT gradient and hessian" begin
    mix = CubicEoS.load(BrusilovskyEoSMixture;
        names=("nitrogen", "methane", "propane", "n-decane"),
    )

    # Two-phase state point
    volume = 1e-6
    nmol = 2000 .* volume .* [0.2463, 0.2208, 0.2208, 0.3121]
    RT = CubicEoS.GAS_CONSTANT_SI * 500

    g1 = similar(nmol, Float64, length(nmol) + 1)
    g1 = CubicEoS.nvtgradient!(g1, mix, nmol, volume, RT)

    g2 = similar(g1)
    h2 = similar(g2, Float64, (size(g2, 1), size(g2, 1)))
    g2, h2 = CubicEoS.nvtgradienthessian!(g2, h2, mix, nmol, volume, RT)

    @test g1 == g2

    h3 = similar(h2)
    h3 = CubicEoS.nvthessian!(h3, mix, nmol, volume, RT)

    @test h2 == h3
end
