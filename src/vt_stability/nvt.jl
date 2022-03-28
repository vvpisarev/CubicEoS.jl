"""
Dimensionless tangent plane distance (TPD) function of Helmholtz free energy.
"""
function helmholtztpd(
    concentration::AbstractVector,
    basestate::VTStabilityBaseState;
    buf::AbstractEoSThermoBuffer=thermo_buffer(basestate.mixture),
)
    gradient = similar(basestate.logcactivity)
    tpd, gradient = helmholtztpdwgradient!(gradient, concentration, basestate; buf=buf)
    return tpd
end

function helmholtztpdwgradient!(
    gradient::AbstractVector,
    concentration::AbstractVector,
    basestate::VTStabilityBaseState;
    buf::AbstractEoSThermoBuffer=thermo_buffer(basestate.mixture),
)
    gradient = helmholtztpdgradient!(gradient, concentration, basestate; buf=buf)
    ptest = pressure(basestate.mixture, concentration, 1, basestate.RT; buf=buf)
    tpd = (-ptest + basestate.pressure) / basestate.RT + dot(concentration, gradient)
    return tpd, gradient
end

"""
[(∇TPD)(c)]ᵢ = (μ'ᵢ - μᵢ) / RT = (ln c'ᵢ + ln γ'ᵢ) - (ln cᵢ + ln γᵢ)
"""
function helmholtztpdgradient!(
    gradient::AbstractVector,
    concentration::AbstractVector,
    basestate::VTStabilityBaseState;
    buf::AbstractEoSThermoBuffer=thermo_buffer(basestate.mixture),
)
    gradient = log_c_activity!(gradient, basestate.mixture, concentration, 1, basestate.RT;
        buf=buf
    )
    gradient .+= log.(concentration)
    gradient .-= basestate.logconcentration
    gradient .-= basestate.logcactivity

    return gradient
end

"""
Hessian of TPD: Hᵢⱼ = δᵢⱼ / c'ᵢ + ∂(ln γᵢ) / ∂c'ⱼ.
"""
function helmholtztpdhessian!(
    hessian::AbstractMatrix,
    concentration::AbstractVector,
    basestate::VTStabilityBaseState;
    buf::AbstractEoSThermoBuffer=thermo_buffer(basestate.mixture),
)
    aux = similar(basestate.logcactivity)  # don't need logs of activity coefficients
    aux, hessian = log_c_activity_wj!(
        aux, hessian, basestate.mixture, concentration, 1, basestate.RT;
        buf=buf
    )
    @inbounds for (ci, hi) in zip(eachindex(concentration), diagind(hessian))
        hessian[hi] += 1 / concentration[ci]
    end
    return hessian
end
