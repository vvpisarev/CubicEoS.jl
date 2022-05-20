function helmholtztpd(
    trialstate::AbstractVTStabilityState,
    basestate::VTStabilityBaseState;
    buf::AbstractEoSThermoBuffer=thermo_buffer(basestate.mixture),
)
    return helmholtztpd(concentration(trialstate), basestate; buf=buf)
end

# This function returns gradient in concentration, not for state variables
# But, I leave it here for future
# function helmholtztpdwgradient!(
#     gradient::AbstractVector,
#     trialstate::AbstractVTStabilityState,
#     basestate::VTStabilityBaseState;
#     buf::AbstractEoSThermoBuffer=thermo_buffer(basestate.mixture),
# )
#     trialconc = concentration(trialstate)
#     tpd, gradient = helmholtztpdwgradient!(gradient, trialconc, basestate; buf=buf)
#     return tpd, gradient
# end

# Works both for Physical and Ideal variables.
"Physical constraint for vt-stability: x + α d ≥ 0."
function physical_constrain_step(
    ::Type{<:AbstractVTStabilityState},
    x::AbstractVector,
    direction::AbstractVector,
)
    αmax = Inf
    αlow = -Inf

    for (i, (xi, di)) in enumerate(zip(x, direction))
        (iszero(di) && xi < 0) && throw(ConstrainStepZeroDirectionError(i, xi))
        αmax = di < 0 ? min(αmax, -xi/di) : αmax
        αlow = di > 0 ? max(αlow, -xi/di) : αlow
    end

    return αlow, αmax
end
