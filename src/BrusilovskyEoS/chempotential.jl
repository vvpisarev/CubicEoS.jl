function CubicEoS.log_c_activity!(
    log_ca::AbstractVector,
    mix::BrusilovskyEoSMixture,
    nmol::AbstractVector,
    volume,
    RT;
    buf::BrusilovskyThermoBuffer=thermo_buffer(mix),
)
    nc = ncomponents(mix)
    ∑mol = sum(nmol)

    comp = components(mix)
    # коэффициенты вычислены для состава в [моль]
    A, B, C, D, aij = eos_parameters(mix, nmol, RT; buf=buf)

    VmB = volume - B
    VpC = volume + C
    VpD = volume + D
    CmD = C - D

    log_VmBbyV = log(VmB/volume)
    ∑molbyVmB = ∑mol / VmB
    AbyRTCmD = A / (RT * CmD)
    log_VpCbyVpD = log(VpC / VpD)

    @inbounds for (i, b, c, d) in zip(1:nc, comp.b, comp.c, comp.d)
        # ∂A_i = 2 ∑_j (N_j * a_ij)
        # прежний код использовал sum( генератор с for in eachindex ) - давал 2 аллокации
        ∂A = 2 * dot(nmol, @view aij[:,i])

        log_ca[i] = -log_VmBbyV + b * ∑molbyVmB
        brackets = (∂A / A - (c - d) / CmD) * log_VpCbyVpD + (c / VpC - d / VpD)
        log_ca[i] -= AbyRTCmD * brackets
    end
    return log_ca
end

function CubicEoS.log_c_activity_wj!(
    log_ca::AbstractVector,
    jacobian::AbstractMatrix,
    mix::BrusilovskyEoSMixture,
    nmol::AbstractVector,
    volume::Real,
    RT::Real;
    buf::BrusilovskyThermoBuffer=thermo_buffer(mix),
)
    A, B, C, D, aij = eos_parameters(mix, nmol, RT; buf=buf)

    aux1 = buf.vec1
    aux2 = buf.vec2
    comp = components(mix)
    nc = ncomponents(mix)
    V = volume
    ntotal = sum(nmol)

    VmB = volume - B
    VpC = volume + C
    VpD = volume + D
    CmD = C - D

    log_VmBbyV = log(VmB/volume)
    AbyRTCmD = A / (RT * CmD)
    log_VpCbyVpD = log(VpC / VpD)
    log2 = 2 * log1p(CmD / (V + D)) / (RT * CmD)

    aux1 .= comp.b ./ VmB

    log_ca .= aux1 .* ntotal .- log_VmBbyV
    jacobian .= aux1 .+ aux1' .+ ntotal .* aux1 .* aux1'
    mul!(aux1, aij, nmol)
    aux2 .= comp.c .- comp.d
    log_ca .-= AbyRTCmD .* log_VpCbyVpD .* (2 .* aux1 ./ A .- aux2 ./ CmD)
    jacobian .-= log2 .* ((A / CmD^2) .* aux2 .* aux2' .+ aij)
    jacobian .+= (log2 / CmD) .* (aux2 .* aux1' .+ aux1 .* aux2')

    aux1 .= 2 .* aux1 ./ (RT * CmD) .- aux2 .* (A / (RT * CmD^2))
    aux2 .= comp.c ./ VpC .- comp.d ./ VpD
    log_ca .-= AbyRTCmD .* aux2
    jacobian .-= aux1 .* aux2'.+ aux2 .* aux1'
    aux1 .= comp.d ./ VpD
    aux2 .= comp.c ./ VpC

    jacobian .-= AbyRTCmD .* (aux1 .* aux1' .- aux2 .* aux2')
    return log_ca, jacobian
end
