"""
    AbstractVTFlashState(x)
    AbstractVTFlashState(concentration, saturation, nmolbase, volumebase)

Abstract type for representation of thermodynamic NVT-state in certain variables `x`.

Second constructor uses `concentration` and `saturation` of a phase and
overall moles and volume of a mixture (`*base`).
"""
abstract type AbstractVTFlashState end

"""
    value(s::AbstractVTFlashState)

Argument for optimization in VT-Flash.
"""
value(s::AbstractVTFlashState) = s.x

"Constructor of a VTFlashState from concentration and saturation."
(::Type{<:AbstractVTFlashState})(concentration, saturation, nmolbase, volumebase) = error("NotImplemented")

"""
    nmolvol!(nmol, state::AbstractVTFlashState, nmolbase, volumebase) -> (nmol, volume)

Destructive option of [`nmolvol`](@ref).
"""
nmolvol!(nmol, s::AbstractVTFlashState, nmolbase, volumebase) = error("NotImplemented")

"""
    nmolvol(state::AbstractVTFlashState, nmolbase, volumebase) -> (nmol, volume)

Converts a `state` variables into moles `nmol` and `volume`
from `state` of a corresponding phase.

See also [`nmolvol!`](@ref).
"""
function nmolvol(s::AbstractVTFlashState, nmolb::V, volumeb::T) where {V, T}
    nmol = similar(nmolb, Float64)  # FIXME: type inference
    return nmolvol!(nmol, s, nmolb, volumeb)
end
