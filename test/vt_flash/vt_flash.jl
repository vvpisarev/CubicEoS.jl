@testset "VT-flash" begin
    mixture = CubicEoS.load(BrusilovskyEoSMixture; names=("methane", "n-pentane"))
    volume = 1e-6

    # Two-phase state point
    nmol = 5000 .* [0.547413, 0.452587] .* volume  # concentration * fraction * volume
    RT = CubicEoS.GAS_CONSTANT_SI * 371

    # Just works test
    @test converged(vt_flash(mixture, nmol, volume, RT))
end
