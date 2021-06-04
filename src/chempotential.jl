#=
Functions to compute chemical potential-related characteristics
=#

"""
    log_activity(mixture, nmol, volume, RT; kwargs...)

Return vector of ln(c_a) - logarithm of activity coefficient
for components of `mixture` at given `nmol`, `volume`, `RT`.
If buffers are provided as keyword arguments, their contents
are modified during the intermediate calculations.

# Arguments

- `mix`: mixture
- `nmol::AbstractVector`: amount of each component (mol)
- `volume`: volume of the mixture (m³)
- `RT`: thermal energy (J/mol)

# Returns

- `AbstractVector`: the logarithms of activity coefficients
of the components at given number of moles, volume and temperature

# Keyword arguments

- `log_a::AbstractVector`: vector for storing the answer. Specify it to avoid allocation
- `ai::AbstractVector`: buffer for intermediate calculations
- `aij::AbstractMatrix`: buffers for intermediate calculations

See (Jirí Mikyska, Abbas Firoozabadi // 10.1002/aic.12387)
"""
function log_activity(
    mix::BrusilovskyEoSMixture{T},
    nmol::AbstractVector,
    volume::Real,
    RT::Real;
    kwargs...
) where {T}

    nc = ncomponents(mix)
    ∑mol = sum(nmol)

    buffers = kwargs.data
    aij = haskey(buffers, :aij) ? buffers[:aij] : Matrix{T}(undef, nc, nc)
    log_a = haskey(buffers, :log_a) ? buffers[:log_a] : Vector{T}(undef, nc)

    # коэффициенты вычислены для состава в [моль]
    A, B, C, D, aij = eos_parameters(mix, nmol, RT; aij = aij, ai = log_a)

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

        log_a[i] = -log_VmBbyV + b * ∑molbyVmB
        brackets = (∂A / A - (c - d)/CmD) * log_VpCbyVpD + (c/VpC - d/VpD)
        log_a[i] -= AbyRTCmD * brackets
    end
    return log_a
end

"""
    log_activity_wj(mixture, nmol, volume, RT; kwargs...)

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

# Returns

- `AbstractVector`: the logarithms of activity coefficients
of the components at given number of moles, volume and temperature

# Keyword arguments

- `log_a::AbstractVector`: vector for storing the answer. Specify it to avoid allocation
- `jacobian::AbstractMatrix`: vector for storing the jacobian. Specify it to avoid allocation
- `aux1::AbstractVector`: buffer for intermediate calculations
- `aux2::AbstractVector`: buffer for intermediate calculations
"""
function log_activity_wj(
    mix::BrusilovskyEoSMixture{T},
    nmol::AbstractVector,
    volume::Real,
    RT::Real;
    kwargs...
) where {T}
    nc = ncomponents(mix)
    buffers = kwargs.data
    jacobian = haskey(buffers, :jacobian) ? buffers[:jacobian] : Matrix{T}(undef, nc, nc)
    log_a = haskey(buffers, :log_a) ? buffers[:log_a] : Vector{T}(undef, nc)
    aux1 = haskey(buffers, :aux1) ? buffers[:aux1] : Vector{T}(undef, nc)
    aux2 = haskey(buffers, :aux2) ? buffers[:aux2] : Vector{T}(undef, nc)

    A, B, C, D, aij = eos_parameters(mix, nmol, RT; aij = jacobian, ai = aux1)

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

    aij = jacobian
    mul!(aux1, aij, nmol)
    @inbounds map!(aux2, 1:nc) do i
        mix[i].c - mix[i].d
    end
    log_a .= (-AbyRTCmD * log_VpCbyVpD) .* (2 .* aux1 ./ A .- aux2 ./ CmD)
    jacobian .= (-2 * log2) .* ((A / CmD^2) .* aux2 .* aux2' .+ aij)
    jacobian .+= (2 * log2 / CmD) .* (aux2 .* aux1' .+ aux1 .* aux2')

    aux1 .= aux1 .* (2 / (RT * CmD)) .- aux2 .* (A / (RT * CmD^2))
    @inbounds map!(aux2, 1:nc) do i
        mix[i].c / VpC - mix[i].d / VpD
    end
    log_a .-= AbyRTCmD .* aux2
    jacobian .-= aux1 .* aux2'.+ aux2 .* aux1'

    @inbounds map!(aux1, 1:nc) do i
        mix[i].b / VmB
    end
    log_a .+= aux1 .* ntotal .- log_VmBbyV
    jacobian .+= aux1 .+ aux1' .+ ntotal .* aux1 .* aux1'

    @inbounds map!(aux1, 1:nc) do i
        mix[i].d / VpD
    end
    @inbounds map!(aux2, 1:nc) do i
        mix[i].c / VpC
    end

    jacobian .-= AbyRTCmD .* (aux1 .* aux1' .- aux2 .* aux2')
    return log_a, jacobian
end
