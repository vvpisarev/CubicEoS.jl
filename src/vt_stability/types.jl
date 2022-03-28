struct VTStabilityBaseState{F<:Real,V<:AbstractVector{<:F},M<:BrusilovskyEoSMixture{<:F}}
    logconcentration::V
    RT::F
    mixture::M
    logcactivity::V
    pressure::F
end

function VTStabilityBaseState(mixture, nmol, volume, RT)
    logconc = log.(nmol ./ volume)
    RT = float(RT)
    logcactivity = log_c_activity(mixture, nmol, volume, RT)
    pres = pressure(mixture, nmol, volume, RT)
    return VTStabilityBaseState(logconc, RT, mixture, logcactivity, pres)
end
