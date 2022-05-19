"""
    AbstractVTFlashState(x)
    AbstractVTFlashState(concentration, saturation, nmolbase, volumebase)

Abstract type for representation of thermodynamic NVT-state in certain variables `x`.

Second constructor uses `concentration` and `saturation` of a phase and
overall moles and volume of a mixture (`*base`).
"""
abstract type AbstractVTFlashState end

struct VTFlashResult{T, S<:AbstractVTFlashState}
    singlephase::Bool
    RT::T
    nmolgas::Vector{T}
    volumegas::T
    nmolliq::Vector{T}
    volumeliq::T
    state::S
    optim::OptimStats
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
    optim = OptimStats(converged, iters, calls)
    return VTFlashResult{T, S}(
        singlephase, RT, nmolgas, volumegas, nmolliq, volumeliq, state, optim
    )
end
