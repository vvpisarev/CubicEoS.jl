"""
    VTStabilityPhysicalState(concentration)

Molar densities (concentrations).
"""
struct VTStabilityPhysicalState{V<:AbstractVector} <: AbstractVTStabilityState
    x::V
end

concentration(state::VTStabilityPhysicalState) = copy(value(state))

function helmholtztpdgradient!(
    gradient,
    trialstate::VTStabilityPhysicalState,
    basestate::VTStabilityBaseState;
    buf::AbstractEoSThermoBuffer=thermo_buffer(basestate.mixture),
)
    trialconc = value(trialstate)
    gradient = helmholtztpdgradient!(gradient, trialconc, basestate; buf=buf)
    return gradient
end

function helmholtztpdhessian!(
    hessian::AbstractMatrix,
    trialstate::VTStabilityPhysicalState,
    basestate::VTStabilityBaseState;
    buf::AbstractEoSThermoBuffer=thermo_buffer(basestate.mixture),
)
    trialconc = value(trialstate)
    return helmholtztpdhessian!(hessian, trialconc, basestate; buf=buf)
end
