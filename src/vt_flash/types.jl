"""
    AbstractVTFlashState(x)
    AbstractVTFlashState(concentration, saturation, nmolb, volumeb)

Abstract type for representation of thermodynamic NVT-state in certain variables `x`.

Second constructor uses concentration and saturation of a phase and
moles and volume of base phase.
"""
abstract type AbstractVTFlashState end

"""
    value(s::AbstractVTFlashState)

Argument for optimization in VT-Flash.
"""
value(s::AbstractVTFlashState) = s.x

"""
    nmolvol(s::AbstractVTFlashState, nmolb, volumeb) -> (nmol, volume)

Moles [mol] and volume [m³] of a phase at `s`tate.
Moles `nmolb` and `volumeb` relate to base state.
"""
nmolvol(s::AbstractVTFlashState, nmolb, volumeb) = error("NotImplemented")

(::Type{<:AbstractVTFlashState})(concentration, saturation, nmolb, volumeb) = error("NotImplemented")

# TODO: deprecate in favor of OptimStats
struct VTFlashOptimStats
    converged::Bool
    iters::Int
    calls::Int
end

struct VTFlashResult{T, S}  # where {S<:AbstractVTFlashState}
    singlephase::Bool
    RT::T
    nmolgas::Vector{T}
    volumegas::T
    nmolliq::Vector{T}
    volumeliq::T
    state::S
    optim::VTFlashOptimStats
end

@inline converged(x::VTFlashResult) = x.optim.converged

function VTFlashResult{T, S}(;
    singlephase,
    converged,
    RT,
    nmolgas,
    volumegas,
    nmolliq,
    volumeliq,
    state=fill(NaN, length(nmolgas)+1),
    iters=-1,
    calls=-1,
) where {T, S}
    optim = VTFlashOptimStats(converged, iters, calls)
    return VTFlashResult{T, S}(
        singlephase, RT, nmolgas, volumegas, nmolliq, volumeliq, state, optim
    )
end
