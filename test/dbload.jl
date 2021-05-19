using Test, CubicEoS, CubicEoSDatabase

"Non-recursively compare field values of `x` and `y`."
function compare_structs(x::T, y::T) where {T}
    for field in fieldnames(T)
        @test getfield(x, field) == getfield(y, field)
    end
end

@testset "src/dbload.jl" begin
    @testset "Default databases" begin
        comp_explicit = load(
            BrusilovskyEoSComponent;
            name="methane",
            component_dbs=(Data.martinez(), Data.brusilovsky_comp())
        )
        comp_default = load(BrusilovskyEoSComponent,
            name="methane"
        )
        compare_structs(comp_explicit, comp_default)

        mix_explicit = load(
            BrusilovskyEoSMixture;
            names = ("methane", "n-pentane"),
            component_dbs=(Data.martinez(), Data.brusilovsky_comp()),
            mix_eos_db = Data.brusilovsky_mix()
        )
        mix_default = load(
            BrusilovskyEoSMixture;
            names = ("methane", "n-pentane"),
        )
        compare_structs(mix_explicit, mix_default)
    end
end
