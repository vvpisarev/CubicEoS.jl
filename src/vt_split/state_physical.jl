"""
    VTSplitPhysicalState(x)
    VTSplitPhysicalState(concentration, saturation, nmolbase, volumebase)

Moles and volume of '-phase variables (a.k.a physical).

**Definition**

```
x = [N'₁, ..., N'ₙ, V']
```

where `N'ᵢ` and `V'`  are moles and volume of `'`-phase, respectively.

See also [`CubicEoS.AbstractVTSplitState`](@ref).
"""
struct VTSplitPhysicalState{V<:AbstractVector} <: AbstractVTSplitState
    x::V
end

function nmolvol!(nmol, s::VTSplitPhysicalState, nmolbase, volumebase)
    x = value(s)
    nmol .= x[1:end-1]
    volume = x[end]
    return nmol, volume
end

function VTSplitPhysicalState{V}(
    concentration::AbstractVector,
    saturation::Real,
    nmolb::AbstractVector,
    volumeb::Real
) where {V}
    x = similar(nmolb, Float64, length(nmolb) + 1)
    volume1 = saturation * volumeb
    @. x[1:end-1] = concentration * volume1
    x[end] = volume1
    return VTSplitPhysicalState{V}(x)
end

@inline VTSplitPhysicalState(c, s, n, v) = VTSplitPhysicalState{Vector{Float64}}(c, s, n, v)

function gradient!(
    grad::AbstractVector,
    state::VTSplitPhysicalState,
    mix,
    nmolb,
    volumeb,
    RT;
    buf::AbstractEoSThermoBuffer=thermo_buffer(mix),
)
    nmol1 = @view state.x[1:end-1]
    volume1 = state.x[end]

    grad = nvtgradient!(grad, mix, nmol1, volume1, RT; buf=buf)

    grad2 = similar(grad)
    nmol2 = similar(nmol1)
    @. nmol2 = nmolb - nmol1
    volume2 = volumeb - volume1
    grad2 = nvtgradient!(grad2, mix, nmol2, volume2, RT; buf=buf)

    grad .-= grad2

    return grad
end

function hessian!(
    hess::AbstractMatrix,
    state::VTSplitPhysicalState,
    mix,
    nmolb,
    volumeb,
    RT;
    buf::AbstractEoSThermoBuffer=thermo_buffer(mix),
)
    nmol1 = @view state.x[1:end-1]
    volume1 = state.x[end]

    hess = nvthessian!(hess, mix, nmol1, volume1, RT; buf=buf)

    hess2 = similar(hess)
    nmol2 = similar(nmol1)
    @. nmol2 = nmolb - nmol1
    volume2 = volumeb - volume1
    hess2 = nvthessian!(hess2, mix, nmol2, volume2, RT; buf=buf)

    hess .+= hess2
    return hess
end

@inline function physical_constrain_step_uplims(
    ::Type{<:VTSplitPhysicalState},
    nmolbase::AbstractVector,
    volumebase::Real,
)
    return [nmolbase; volumebase]
end
