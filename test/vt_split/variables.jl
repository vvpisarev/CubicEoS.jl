@testset "Different variables" begin
    mix = CubicEoS.load(BrusilovskyEoSMixture;
        names=("nitrogen", "methane", "propane", "n-decane"),
    )

    # Two-phase state point
    volume = 1e-6
    nmol = 2000 .* volume .* [0.2463, 0.2208, 0.2208, 0.3121]
    RT = CubicEoS.GAS_CONSTANT_SI * 500

    # Hard-coded concentration found from vt-stability for that mixture and negative tpd.
    trialconc = [510, 432, 372, 332]

    for vars in CubicEoS.vtsplitvariables()
        @testset "Flash converged for $vars" begin
            @test converged(
                vt_split(mix, nmol, volume, RT, trialconc, vars;
                    eos_constrain_step=CubicEoS.BrusilovskyEoS.eos_vt_split_constrain_step(vars, mix, nmol, volume),
                )
            )
        end
    end

    @test_broken vt_split_newton(mix, nmol, volume, RT, trialconc, vars; tol=1e-8)

    # @testset "Equality of approaches" begin
    #     resultbfgs = map(StateTypes) do ST
    #         vt_split(mix, nmol, volume, RT, ST; tol=1e-8)
    #     end
    #     resultnewton = map(StateTypes) do ST
    #         vt_split_newton(mix, nmol, volume, RT, ST; tol=1e-8)
    #     end
    #     results = tuple(resultbfgs..., resultnewton...)

    #     function compare(x, y)
    #         @test x.nmolgas ≈ y.nmolgas rtol=1e-5
    #         @test x.volumegas ≈ y.volumegas rtol=1e-5
    #         @test x.nmolliq ≈ y.nmolliq rtol=1e-5
    #         @test x.volumeliq ≈ y.volumeliq rtol=1e-5
    #     end
    #     for (r1, r2) in zip(results[1:end-1], results[2:end])
    #         compare(r1, r2)
    #     end
    #     for (rbfgs, rnewt) in zip(resultbfgs, resultnewton)
    #         @test CubicEoS.value(rbfgs.state) ≈ CubicEoS.value(rnewt.state) rtol=1e-5
    #     end
    # end
end
