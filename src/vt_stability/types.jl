# BrusilovskyEoSMixture -> AbstractEoSMixture
struct VTStabilityBaseState{F<:Real,V<:AbstractVector{<:F},M<:AbstractEoSMixture}
    logconcentration::V
    RT::F
    mixture::M
    logcactivity::V
    pressure::F
end

function VTStabilityBaseState(mixture, nmol, volume, RT; buf=thermo_buffer(mixture))
    logconc = log.(nmol ./ volume)
    RT = float(RT)
    logcactivity = log_c_activity(mixture, nmol, volume, RT; buf=buf)
    pres = pressure(mixture, nmol, volume, RT; buf=buf)
    return VTStabilityBaseState(logconc, RT, mixture, logcactivity, pres)
end

abstract type AbstractVTStabilityState end

value(state::AbstractVTStabilityState) = state.x
concentration(state::AbstractVTStabilityState) = error("Not implemented")
concentration!(conc, ::Type{<:AbstractVTStabilityState}, stateval) = error("Not implemented")
fromconcentration(::Type{<:AbstractVTStabilityState}, x) = error("Not implemented")

struct VTStabilityResult{T<:Real,S<:AbstractVTStabilityState}
    issuccess::Bool
    isstable::Bool
    energy_density::T
    concentration::Vector{T}
    state::S
    optim::OptimStats
end
