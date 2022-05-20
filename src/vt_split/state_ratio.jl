"""
    VTSplitRatioState(x)
    VTSplitRatioState(concentration, saturation, nmolbase, volumebase)

Fractional (ratio) variables for flash.
Moles and volume of '-phase variables (a.k.a physical).

**Definition**

```
x = [N'₁/N₁, ..., N'ₙ/Nₙ, V'/V]
```

where

- `N'ᵢ` and `V'`  are moles and volume of `'`-phase, respectively;
- `Nᵢ` and `V` are moles and volume of mixture (`nmolbase`, `volumebase`), respectively.

See also [`CubicEoS.AbstractVTSplitState`](@ref).
"""
struct VTSplitRatioState{V<:AbstractVector} <: AbstractVTSplitState
    x::V
end

function nmolvol!(nmol, s::VTSplitRatioState, nmolbase, volumebase)
    x = value(s)
    @. nmol = nmolbase .* x[1:end-1]
    volume = volumebase * x[end]
    return nmol, volume
end

function VTSplitRatioState{V}(
    concentration::AbstractVector,
    saturation::Real,
    nmolb::AbstractVector,
    volumeb::Real
) where {V}
    x = similar(nmolb, Float64, length(nmolb) + 1)
    @. x[1:end-1] = volumeb * saturation * concentration / nmolb
    x[end] = saturation
    return VTSplitRatioState{V}(x)
end

@inline VTSplitRatioState(c, s, n, v) = VTSplitRatioState{Vector{Float64}}(c, s, n, v)

function gradient!(
    grad::AbstractVector,
    state::VTSplitRatioState,
    mix::AbstractEoSMixture,
    nmolb::AbstractVector,
    volumeb::Real,
    RT::Real;
    buf::AbstractEoSThermoBuffer=thermo_buffer(mix),
)
    nmol1, vol1 = nmolvol(state, nmolb, volumeb)
    grad = nvtgradient!(grad, mix, nmol1, vol1, RT; buf=buf)

    nmol2 = nmolb - nmol1  # alloc
    vol2 = volumeb - vol1
    grad2 = similar(grad)  # alloc
    grad2 = nvtgradient!(grad2, mix, nmol2, vol2, RT; buf=buf)

    grad .-= grad2
    grad[1:end-1] .*= nmolb
    grad[end] *= volumeb

    return grad
end

function hessian!(
    hess::AbstractMatrix,
    state::VTSplitRatioState,
    mix::AbstractEoSMixture,
    nmolb::AbstractVector,
    volumeb::Real,
    RT::Real;
    buf::AbstractEoSThermoBuffer=thermo_buffer(mix),
)
    nmol1, vol1 = nmolvol(state, nmolb, volumeb)
    hess = nvthessian!(hess, mix, nmol1, vol1, RT; buf=buf)

    nmol2 = nmolb - nmol1  # alloc
    vol2 = volumeb - vol1
    hess2 = similar(hess)  # alloc
    hess2 = nvthessian!(hess2, mix, nmol2, vol2, RT; buf=buf)

    hess .+= hess2

    # Also, kron and gemm! is applicable
    scale = [nmolb; volumeb] * [nmolb; volumeb]'  # alloc
    hess .*= scale

    return hess
end

@inline function physical_constrain_step_uplims(
    ::Type{<:VTSplitRatioState},
    nmolbase::AbstractVector,
    volumebase::Real=NaN,
)
    return ones(size(nmolbase, 1) + 1)
end
