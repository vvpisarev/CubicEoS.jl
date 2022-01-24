@testset "NVT gradient and hessian" begin
    mix = CubicEoS.load(BrusilovskyEoSMixture; names=("methane", "n-pentane"))
    volume = 1e-6

    # Two-phase state point
    nmol = 5000 .* [0.547413, 0.452587] .* volume  # concentration * fraction * volume
    RT = CubicEoS.GAS_CONSTANT_SI * 371

    g1 = similar(nmol, Float64, length(nmol) + 1)
    g1 = CubicEoS.nvtgradient!(g1, mix, nmol, volume, RT)

    g2 = similar(g1)
    h2 = similar(g2, Float64, (size(g2, 1), size(g2, 1)))
    g2, h2 = CubicEoS.nvtgradienthessian!(g2, h2, mix, nmol, volume, RT)

    @test g1 == g2
end
