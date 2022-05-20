"""
    AbstractVTSplitState(x)
    AbstractVTSplitState(concentration, saturation, nmolbase, volumebase)

Abstract type for representation of thermodynamic NVT-state in certain variables `x`.

Second constructor uses `concentration` and `saturation` of a phase and
overall moles and volume of a mixture (`*base`).
"""
abstract type AbstractVTSplitState end

struct VTSplitResult{T, S<:AbstractVTSplitState}
    singlephase::Bool
    RT::T
    nmolgas::Vector{T}
    volumegas::T
    nmolliq::Vector{T}
    volumeliq::T
    state::S
    optim::OptimStats
end

converged(x::VTSplitResult) = converged(x.optim)
iters(x::VTSplitResult) = iters(x.optim)
calls(x::VTSplitResult) = calls(x.optim)

function VTSplitResult{T, S}(;
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
    return VTSplitResult{T, S}(
        singlephase, RT, nmolgas, volumegas, nmolliq, volumeliq, state, optim
    )
end
