abstract type AbstractEoSComponent end
abstract type AbstractEoSMixture{T} end
abstract type AbstractEoSThermoBuffer end

# TODO: move to interface.jl
"""
    thermo_buffer(mix)

Create a buffer for intermediate calculations of mixture thermodynamic properties.

See also: [`pressure`](@ref), [`log_c_activity`](@ref), [`log_c_activity!`](@ref),
[`log_c_activity_wj`](@ref), [`log_c_activity_wj!`](@ref)
"""
thermo_buffer(::AbstractEoSMixture) = error("NotImplemented")
