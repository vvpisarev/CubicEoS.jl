@testset "vt_flash" begin
    mix = CubicEoS.load(BrusilovskyEoSMixture;
        names=("nitrogen", "methane", "propane", "n-decane"),
    )

    # Two-phase state point
    volume = 1e-6
    nmol = 2000 .* volume .* [0.2463, 0.2208, 0.2208, 0.3121]
    RT = CubicEoS.GAS_CONSTANT_SI * 500

    StateTypes = (
        CubicEoS.PhysicalState,
        CubicEoS.RatioState,
        CubicEoS.IdealIdentityState,
    )

    # Just works tests
    @testset "Converged" begin
        @test converged(vt_flash(mix, nmol, volume, RT))
        for ST in StateTypes
            @test converged(vt_flash(mix, nmol, volume, RT, ST))
        end
    end

    @testset "Equality of approaches" begin
        resultbfgs = map(StateTypes) do ST
            vt_flash(mix, nmol, volume, RT, ST; gtol=1e-5/RT)
        end
        resultnewton = map(StateTypes) do ST
            vt_flash_newton(mix, nmol, volume, RT, ST; gtol=1e-5/RT)
        end
        results = tuple(resultbfgs..., resultnewton...)

        function compare(x, y)
            @test x.nmolgas ≈ y.nmolgas rtol=1e-5
            @test x.volumegas ≈ y.volumegas rtol=1e-5
            @test x.nmolliq ≈ y.nmolliq rtol=1e-5
            @test x.volumeliq ≈ y.volumeliq rtol=1e-5
        end
        for (r1, r2) in zip(results[1:end-1], results[2:end])
            compare(r1, r2)
        end
        for (rbfgs, rnewt) in zip(resultbfgs, resultnewton)
            @test CubicEoS.value(rbfgs.state) ≈ CubicEoS.value(rnewt.state) rtol=1e-5
        end
    end
end
