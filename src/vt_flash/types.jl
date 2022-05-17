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
