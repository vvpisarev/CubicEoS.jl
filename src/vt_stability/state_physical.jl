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

function __constrain_step(
    ::Type{VTStabilityPhysicalState},
    trialx::AbstractVector,
    direction::AbstractVector,
    covolumes::AbstractVector,
)
    αmax = Inf
    αlow = -Inf

    # non-negativness
    for (i, (ci, di)) in enumerate(zip(trialx, direction))
        (iszero(di) && ci < 0) && throw(ConstrainStepZeroDirectionError(i, ci))
        αmax = di < 0 ? min(αmax, -ci/di) : αmax
        αlow = di > 0 ? max(αlow, -ci/di) : αlow
    end

    # size of phase
    dirdotcov = dot(direction, covolumes)
    trialxdotcov = dot(trialx, covolumes)
    # case `iszero(dirdotcov)` is captured too
    αmax = dirdotcov > 0 ? min(αmax, (1 - trialxdotcov)/dirdotcov) : αmax
    αlow = dirdotcov < 0 ? max(αlow, (1 - trialxdotcov)/dirdotcov) : αlow

    (αmax ≤ 0 || αmax ≤ αlow) && throw(ConstrainStepLowerBoundError(trialx, direction))

    return αmax
end
