@testset "dbload.jl" begin
    @testset "Default databases" begin
        comp_explicit = CubicEoS.load(BrusilovskyEoSComponent;
            name="methane",
            component_dbs=(Data.martinez(), Data.brusilovsky_comp())
        )
        comp_default = CubicEoS.load(BrusilovskyEoSComponent,
            name="methane",
        )
        compare_structs(comp_explicit, comp_default)

        mix_explicit = CubicEoS.load(BrusilovskyEoSMixture;
            names = ("methane", "n-pentane"),
            component_dbs=(Data.martinez(), Data.brusilovsky_comp()),
            mix_eos_db = Data.brusilovsky_mix(),
        )
        mix_default = CubicEoS.load(BrusilovskyEoSMixture;
            names = ("methane", "n-pentane"),
        )
        compare_structs(mix_explicit, mix_default)
    end
end
