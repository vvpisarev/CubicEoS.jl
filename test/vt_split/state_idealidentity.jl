@testset "VTSplitIdealIdentityState variables" begin
    nmol = [1.5, 2.5, 4.3]
    volume = 0.3e-6

    nmolb = [2.1, 3.0, 4.5]
    volumeb = volume + 0.7e-6

    let concentration = nmol / volume, saturation = volume / volumeb
        state = CubicEoS.VTSplitIdealIdentityState(concentration, saturation, nmolb, volumeb)
        nmolt, volt = CubicEoS.nmolvol(state, nmolb, volumeb)

        @test nmolt ≈ nmol
        @test volume ≈ volt
    end
end
