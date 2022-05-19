@testset "VT-stability" begin
    @testset "Minimal Working Example" begin
        mix = CubicEoS.load(BrusilovskyEoSMixture;
            names=("nitrogen", "methane", "propane", "n-decane"),
        )

        # Two-phase state point
        volume = 1e-6
        nmol = 2000 .* volume .* [0.2463, 0.2208, 0.2208, 0.3121]
        RT = CubicEoS.GAS_CONSTANT_SI * 500

        issuccess, isstable, results = vt_stability(mix, nmol, volume, RT)

        @test issuccess
        @test !isstable
    end
end
