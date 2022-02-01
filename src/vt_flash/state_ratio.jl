#=
Flash based on fractional (ratio) variables

        [ N'1       N'n   V' ]
state = [ ---, ..., ---, --- ]
        [  N1        Nn   V  ]
=#

struct RatioState{V<:AbstractVector} <: AbstractVTFlashState
    x::V
end

function nmolvol(s::RatioState, nmolb::V, volumeb::T) where {V, T}
    x = value(s)
    nmol1 = nmolb .* x[1:end-1]
    volume1 = volumeb * x[end]
    return (nmol1, volume1)
end

function RatioState(
    concentration::AbstractVector,
    saturation::Real,
    nmolb::AbstractVector,
    volumeb::Real
)
    x = similar(nmolb, Float64, length(nmolb) + 1)
    @. x[1:end-1] = volumeb * saturation * concentration / nmolb
    x[end] = saturation
    return RatioState(x)
end
