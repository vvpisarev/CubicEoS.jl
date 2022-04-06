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
