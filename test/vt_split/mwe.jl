@testset "Minimal Working Example" begin
    mix = CubicEoS.load(BrusilovskyEoSMixture;
        names=("nitrogen", "methane", "propane", "n-decane"),
    )

    # Two-phase state point, known a priori
    volume = 1e-6
    nmol = 2000 .* volume .* [0.2463, 0.2208, 0.2208, 0.3121]
    RT = CubicEoS.GAS_CONSTANT_SI * 500

    vars = CubicEoS.VTSplitIdealIdentityState
    eos_constrain = CubicEoS.BrusilovskyEoS.eos_vt_split_constrain_step(vars, mix, nmol, volume)

    _, _, stability_tries = vt_stability(mix, nmol, volume, RT)
    trialconc = CubicEoS.concentrationwithlowesttpd(stability_tries)

    @test let
        result = vt_split(mix, nmol, volume, RT, trialconc, vars;
                eos_constrain_step=eos_constrain,
        )
        converged(result)
    end
end
