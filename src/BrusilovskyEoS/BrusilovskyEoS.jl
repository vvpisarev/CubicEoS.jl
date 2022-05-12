module BrusilovskyEoS
    export BrusilovskyEoSComponent, BrusilovskyEoSMixture

    import ..CubicEoS
    using CubicEoS.CubicEoSDatabase

    # Useful names for dispatch, used as shortcuts for CubicEoS.foo -> foo
    using CubicEoS: AbstractEoSComponent
    using CubicEoS: AbstractEoSMixture
    using CubicEoS: AbstractEoSThermoBuffer, thermo_buffer
    using CubicEoS: ncomponents, components

    using LinearAlgebra
    using StructArrays

    include("types.jl")
    include("interface.jl")
    include("dbload.jl")
    include("basic_thermo.jl")
    include("chempotential.jl")
    include("vt_stability.jl")
end
