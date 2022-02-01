#=
Functions to compute chemical potential-related characteristics
=#

"""
    log_c_activity(mixture, nmol, volume, RT[; buf])

Return vector of ln(c_a) - logarithm of activity coefficients
for components of `mixture` at given `nmol`, `volume`, `RT`.
If buffers are provided as keyword arguments, their contents
are modified during the intermediate calculations.

# Arguments

- `mix`: mixture
- `nmol::AbstractVector`: amount of each component (mol)
- `volume`: volume of the mixture (m³)
- `RT`: thermal energy (J/mol)

# Keywords

- `buf::Union{NamedTuple, AbstractDict, BrusilovskyThermoBuffer}`: buffers for intermediate
    calculations (see [`pressure`](@ref))

# Returns

- `AbstractVector`: the logarithms of activity coefficients of the components at given
    number of moles, volume and temperature

See (Jirí Mikyska, Abbas Firoozabadi // 10.1002/aic.12387)

See also: [`log_c_activity!`](@ref), [`log_c_activity_wj`](@ref), [`log_c_activity_wj!`](@ref)
"""
function log_c_activity(
    mix::BrusilovskyEoSMixture{Tmix},
    nmol::AbstractVector{Tmol},
    volume::Tvol,
    RT::Ttemp;
    buf = NamedTuple(),
) where {Tmix, Tmol, Tvol, Ttemp}

    T = promote_type(Tmix, Tmol, Tvol, Ttemp)
    nc = ncomponents(mix)
    log_ca = Vector{T}(undef, nc)

    return log_c_activity!(log_ca, mix, nmol, volume, RT; buf = buf)
end

"""
    log_c_activity!(log_ca, mixture, nmol, volume, RT[; buf])

Return vector of ln(c_a) - logarithm of activity coefficients
for components of `mixture` at given `nmol`, `volume`, `RT`.
If buffers are provided as keyword arguments, their contents
are modified during the intermediate calculations. The answer is stored
in `log_ca`.

# Arguments

- `log_ca::AbstractVector`: buffer to store the result
- `mix::BrusilovskyEoSMixture`: mixture
- `nmol::AbstractVector`: amount of each component (mol)
- `volume`: volume of the mixture (m³)
- `RT`: thermal energy (J/mol)

# Keywords

- `buf::Union{BrusilovskyThermoBuffer, NamedTuple, AbstractDict}`: buffers for intermediate
    calculations (see [`pressure`](@ref))

# Returns

- `AbstractVector`: the logarithms of activity coefficients of the components at given
    number of moles, volume and temperature

See (Jirí Mikyska, Abbas Firoozabadi // 10.1002/aic.12387)

See also: [`log_c_activity`](@ref), [`log_c_activity_wj`](@ref), [`log_c_activity_wj!`](@ref)
"""
function log_c_activity!(
    log_ca::AbstractVector,
    mix::BrusilovskyEoSMixture,
    nmol::AbstractVector,
    volume,
    RT;
    buf = NamedTuple(),
)
    nc = ncomponents(mix)
    ∑mol = sum(nmol)

    # коэффициенты вычислены для состава в [моль]
    A, B, C, D, aij = eos_parameters(mix, nmol, RT; buf = buf)

    VmB = volume - B
    VpC = volume + C
    VpD = volume + D
    CmD = C - D

    log_VmBbyV = log(VmB/volume)
    ∑molbyVmB = ∑mol / VmB
    AbyRTCmD = A / (RT * CmD)
    log_VpCbyVpD = log(VpC / VpD)

    @inbounds for i in 1:nc
        subst = mix.components[i]
        b = subst.b
        c = subst.c
        d = subst.d

        # ∂A_i = 2 ∑_j (N_j * a_ij)
        # прежний код использовал sum( генератор с for in eachindex ) - давал 2 аллокации
        ∂A = 2 * dot(nmol, @view aij[:,i])

        log_ca[i] = -log_VmBbyV + b * ∑molbyVmB
        brackets = (∂A / A - (c - d) / CmD) * log_VpCbyVpD + (c / VpC - d / VpD)
        log_ca[i] -= AbyRTCmD * brackets
    end
    return log_ca
end

"""
    log_c_activity_wj(mixture, nmol, volume, RT[; buf])

Return vector of ln(c_a) - logarithm of activity coefficient
for components of `mixture` at given `nmol`, `volume`, `RT` -
and the jacobian ∂ln(c_a[i]) / ∂n[j].
If buffers are provided as keyword arguments, their contents
are modified during the intermediate calculations.

# Arguments

- `mix`: mixture
- `nmol::AbstractVector`: amount of each component (mol)
- `volume`: volume of the mixture (m³)
- `RT`: thermal energy (J/mol)

# Keywords

- `buf::BrusilovskyThermoBuffer`: buffers for intermediate calculations
    (see [`thermo_buffer`](@ref))

# Returns

- `Tuple{AbstractVector,AbstractMatrix}`: the logarithms of activity coefficients
of the components at given number of moles, volume and temperature, and the jacobian
matrix ∂ln(c_a[i]) / ∂n[j].

See also: [`log_c_activity`](@ref), [`log_c_activity!`](@ref), [`log_c_activity_wj!`](@ref)
"""
function log_c_activity_wj(
    mix::BrusilovskyEoSMixture{T},
    nmol::AbstractVector,
    volume::Real,
    RT::Real;
    buf::BrusilovskyThermoBuffer = thermo_buffer(mix),
) where {T}
    nc = ncomponents(mix)
    log_ca = Vector{T}(undef, nc)
    jacobian = similar(log_ca, (nc, nc))
    return log_c_activity_wj!(log_ca, jacobian, mix, nmol, volume, RT; buf = buf)
end

"""
    log_activity_wj!(log_ca, jacobian, mixture, nmol, volume, RT[; buf])

Return vector of ln(c_a) - logarithm of activity coefficient
for components of `mixture` at given `nmol`, `volume`, `RT` -
and the jacobian ∂ln(c_a[i]) / ∂n[j]. The first two arguments get overwritten
by the result.
If buffers are provided as keyword arguments, their contents
are modified during the intermediate calculations.

# Arguments

- `log_ca::AbstractVector`: a vector to store the activity coefficients
- `jacobian::AbstractMatrix`: a matrix to store ∂ln(c_a[i]) / ∂n[j]
- `mix`: mixture
- `nmol::AbstractVector`: amount of each component (mol)
- `volume`: volume of the mixture (m³)
- `RT`: thermal energy (J/mol)

# Keywords

- `buf::BrusilovskyThermoBuffer`: buffers for intermediate calculations
    (see [`thermo_buffer`](@ref))

# Returns

- `Tuple{AbstractVector,AbstractMatrix}`: the logarithms of activity coefficients
of the components at given number of moles, volume and temperature, and the jacobian
matrix ∂ln(c_a[i]) / ∂n[j]. The values are aliases for the first two arguments.

See also: [`log_c_activity`](@ref), [`log_c_activity!`](@ref), [`log_c_activity_wj`](@ref)
"""
function log_c_activity_wj!(
    log_ca::AbstractVector,
    jacobian::AbstractMatrix,
    mix::BrusilovskyEoSMixture,
    nmol::AbstractVector,
    volume::Real,
    RT::Real;
    buf::BrusilovskyThermoBuffer = thermo_buffer(mix),
)
    return __log_c_activity_wj_impl__(log_ca, jacobian, mix, nmol, volume, RT, buf)
end

function __log_c_activity_wj_impl__(
    log_ca::AbstractVector,
    jacobian::AbstractMatrix,
    mix::BrusilovskyEoSMixture,
    nmol::AbstractVector,
    volume,
    RT,
    buf::BrusilovskyThermoBuffer,
)
    A, B, C, D, aij = eos_parameters(mix, nmol, RT; buf = buf)

    aux1 = buf.vec1
    aux2 = buf.vec2
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
    log2 = log1p(CmD / (V + D)) / (RT * CmD)

    @inbounds map!(aux1, 1:nc) do i
        mix[i].b / VmB
    end
    log_ca .= aux1 .* ntotal .- log_VmBbyV
    jacobian .= aux1 .+ aux1' .+ ntotal .* aux1 .* aux1'
    mul!(aux1, aij, nmol)
    @inbounds map!(aux2, 1:nc) do i
        mix[i].c - mix[i].d
    end
    log_ca .-= AbyRTCmD .* log_VpCbyVpD .* (2 .* aux1 ./ A .- aux2 ./ CmD)
    jacobian .-= (2 * log2) .* ((A / CmD^2) .* aux2 .* aux2' .+ aij)
    jacobian .+= (2 * log2 / CmD) .* (aux2 .* aux1' .+ aux1 .* aux2')

    aux1 .= 2 .* aux1 ./ (RT * CmD) .- aux2 .* (A / (RT * CmD^2))
    @inbounds map!(aux2, 1:nc) do i
        mix[i].c / VpC - mix[i].d / VpD
    end
    log_ca .-= AbyRTCmD .* aux2
    jacobian .-= aux1 .* aux2'.+ aux2 .* aux1'
    @inbounds map!(aux1, 1:nc) do i
        mix[i].d / VpD
    end
    @inbounds map!(aux2, 1:nc) do i
        mix[i].c / VpC
    end

    jacobian .-= AbyRTCmD .* (aux1 .* aux1' .- aux2 .* aux2')
    return log_ca, jacobian
end
