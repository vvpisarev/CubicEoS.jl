#=
    IdealIdentityState

Variables that make ideal part of hessian over moles to identity matrix.

state = [ 2 √Nᵢ arcsin(√N'ᵢ / √Nᵢ); √N V'/V ], where N = Σᵢ Nᵢ (total moles).
         |-----------------------| |-------|
                    α'ᵢ               α'V
=#

struct IdealIdentityState{V<:AbstractVector} <: AbstractVTFlashState
    x::V
end

function nmolvol(s::IdealIdentityState, nmolb::V, volumeb::T) where {V, T}
    x = value(s)
    αnmol = @view x[1:end-1]
    αvol = x[end]

    nmol = @. nmolb * (sin(αnmol/(2 * sqrt(nmolb))))^2

    volume = (αvol * volumeb) / sqrt(sum(nmolb))
    return (nmol, volume)
end

function IdealIdentityState(
    concentration::AbstractVector,
    saturation::Real,
    nmolb::AbstractVector,
    volumeb::Real,
)
    vol = saturation * volumeb
    nmol = concentration * vol

    x = Vector{Float64}(undef, length(nmolb) + 1)

    @. x[1:end-1] = 2 * sqrt(nmolb) * asin(sqrt(nmol/nmolb))
    x[end] = sqrt(sum(nmolb)) * saturation

    return IdealIdentityState(x)
end
