"""
    VTStabilityIdealIdentityState(x)

The variables are `2 * sqrt(concentration)`.
This choice makes ideal part hessian identity matrix.

**Reference**
1. Nichita D.V., Fluid Ph. Equil., 2017, 10.1016/j.fluid.2017.05.022
2. Michelsen M.L., Fluid Ph. Equil., 1982, 10.1016/0378-3812(82)85001-2
"""
struct VTStabilityIdealIdentityState{V<:AbstractVector} <: AbstractVTStabilityState
    x::V
end

function concentration!(c, ::Type{<:VTStabilityIdealIdentityState}, stval)
    @. c = stval^2 / 4
    return c
end

concentration(state::VTStabilityIdealIdentityState) = value(state).^2 ./ 4

function fromconcentration(::Type{<:VTStabilityIdealIdentityState}, x)
    vars = 2 .* sqrt.(x)
    return VTStabilityIdealIdentityState(vars)
end

"""
Transition from concentration-variables is [(∇TPD)(α)]ᵢ = [(∇TPD)(c)]ᵢ * √cᵢ
"""
function helmholtztpdgradient!(
    gradient::AbstractVector,
    trialstate::VTStabilityIdealIdentityState,
    basestate::VTStabilityBaseState;
    buf::AbstractEoSThermoBuffer=thermo_buffer(basestate.mixture),
)
    trialconc = concentration(trialstate)
    gradient = helmholtztpdgradient!(gradient, trialconc, basestate; buf=buf)
    gradient .*= sqrt.(trialconc)
    return gradient
end

"""
Transition from concentration-variables is

[(∇²TPD)(α)]ᵢⱼ = [(∇²TPD)(c)]ᵢⱼ √cᵢ √cⱼ + δᵢⱼ [(∇TPD)(c)]ᵢ / 2
"""
function helmholtztpdhessian!(
    hessian::AbstractMatrix,
    trialstate::VTStabilityIdealIdentityState,
    basestate::VTStabilityBaseState;
    buf::AbstractEoSThermoBuffer=thermo_buffer(basestate.mixture),
)
    trialconc = concentration(trialstate)
    gradient = similar(trialconc)

    # Concentration-dependant gradient and Hessian
    gradient, hessian = helmholtztpdgradienthessian!(gradient, hessian, trialconc, basestate;
        buf=buf
    )

    # Hessian scaling
    @. hessian *= sqrt(trialconc * trialconc')
    @inbounds for (gi, hi) in zip(eachindex(gradient), diagind(hessian))
        hessian[hi] += gradient[gi] / 2
    end

    return hessian
end
