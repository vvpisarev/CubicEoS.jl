#=
Functions to compute chemical potential-related characteristics
=#

struct ChemPotentialMixture{M, B}<:AbstractEoSMixture
    mixture_eos::M
    buffers::B
end

function ChemPotentialMixture(mix::BrusilovskyEoSMixture{T}) where {T}
    nc = ncomponents(mix)
    buffers = (
        aij = zeros(T, (nc, nc)),
        jacobian = zeros(T, (nc, nc)),
        log_a = zeros(T, nc),
        aux1 = zeros(T, nc),
        aux2 = zeros(T, nc)
    )
    return ChemPotentialMixture(mix, buffers)
end

function pressure!(
    mix::ChemPotentialMixture{<:BrusilovskyEoSMixture},
    nmol::AbstractVector,
    volume,
    RT
)
    buf = mix.buffers
    ai = buf.aux1
    aij = buf.aij
    return pressure!(aij, ai, mix.mixture_eos, nmol, volume, RT)
end

"""
    log_activity!([log_a_, ai_, aij_,] mixture, nmol, volume, RT)

Returns vector of ln(c_a) - logarithm of activity coefficient
for components of `mixture` at given

- amount of each component `N` (mol)
- volume `V` (m³)
- thermal energy `RT`

# Optional arguments

- `log_a_::AbstractVector` - vector for storing answer. Specify it to avoid allocation
- `ai::AbstractVector`, `aij::AbstractMatrix` - buffer for intermediate calculations

See (Jirí Mikyska, Abbas Firoozabadi // 10.1002/aic.12387)
"""
function log_activity!(
    log_a_::AbstractVector,
    ai_::AbstractVector,
    aij_::AbstractMatrix,
    mix::BrusilovskyEoSMixture{T},
    nmol::AbstractVector,
    volume::Real,
    RT::Real
) where {T}

    nc = ncomponents(mix)
    ∑mol = sum(nmol)

    # коэффициенты вычислены для состава в [моль]
    A, B, C, D, aij = eos_parameters!(aij_, ai_, mix, nmol, RT)

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

        log_a_[i] = -log_VmBbyV + b * ∑molbyVmB
        brackets = (∂A / A - (c - d)/CmD) * log_VpCbyVpD + (c/VpC - d/VpD)
        log_a_[i] -= AbyRTCmD * brackets
    end
    return log_a_
end

function log_activity!(
    mix::ChemPotentialMixture{<:BrusilovskyEoSMixture},
    nmol::AbstractVector,
    volume::Real,
    RT::Real
)

    aij = mix.buffers.aij
    ai = mix.buffers.aux1
    log_a = mix.buffers.log_a
    return log_activity!(log_a, ai, aij, mix.mixture_eos, nmol, volume, RT)
end

function log_activity_wj!(
    log_a_::AbstractVector,
    jacobian::AbstractMatrix,
    aij_::AbstractMatrix,
    aux1::AbstractVector,
    aux2::AbstractVector,
    mix::BrusilovskyEoSMixture{T},
    nmol::AbstractVector,
    volume::Real,
    RT::Real
) where {T}
    A, B, C, D, aij = eos_parameters!(aij_, aux1, mix, nmol, RT)

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
    log_a_ .= aux1 .* ntotal .- log_VmBbyV
    jacobian .= aux1 .+ aux1' .+ ntotal .* aux1 .* aux1'
    mul!(aux1, aij, nmol)
    @inbounds map!(aux2, 1:nc) do i
        mix[i].c - mix[i].d
    end
    log_a_ .-= AbyRTCmD .* log_VpCbyVpD .* (2 .* aux1 ./ A .- aux2 ./ CmD)
    jacobian .-= (2 * log2) .* ((A / CmD^2) .* aux2 .* aux2' .+ aij)
    jacobian .+= (2 * log2 / CmD) .* (aux2 .* aux1' .+ aux1 .* aux2')

    aux1 .= 2 .* aux1 ./ (RT * CmD) .- aux2 .* (A / (RT * CmD^2))
    @inbounds map!(aux2, 1:nc) do i
        mix[i].c / VpC - mix[i].d / VpD
    end
    log_a_ .-= AbyRTCmD .* aux2
    jacobian .-= aux1 .* aux2'.+ aux2 .* aux1'
    @inbounds map!(aux1, 1:nc) do i
        mix[i].d / VpD
    end
    @inbounds map!(aux2, 1:nc) do i
        mix[i].c / VpC
    end

    jacobian .-= AbyRTCmD .* (aux1 .* aux1' .- aux2 .* aux2')
    return log_a_, jacobian
end

function log_activity_wj!(
    mix::ChemPotentialMixture{<:BrusilovskyEoSMixture},
    nmol::AbstractVector,
    volume::Real,
    RT::Real
)
    aij = mix.buffers.aij
    jacobian = mix.buffers.jacobian
    aux1 = mix.buffers.aux1
    aux2 = mix.buffers.aux2
    log_a = mix.buffers.log_a
    return log_activity_wj!(
        log_a,
        jacobian,
        aij,
        aux1,
        aux2,
        mix.mixture_eos,
        nmol,
        volume,
        RT,
    )
end

function log_activity(
    mix::BrusilovskyEoSMixture{T},
    nmol::AbstractVector,
    volume,
    RT
) where {T}
    nc = length(mix)
    aij = zeros(T, nc, nc)
    ai = zeros(T, nc)
    log_a = zeros(T, nc)
    return log_activity!(log_a, ai, aij, mix, nmol, volume, RT)
end

function log_activity_wj(
    mix::BrusilovskyEoSMixture{T},
    nmol::AbstractVector,
    volume,
    RT
) where {T}
    nc = length(mix)
    aij = zeros(T, nc, nc)
    jacobian = zeros(T, nc, nc)
    aux1, aux2 = zeros(T, nc), zeros(T, nc)
    log_a = zeros(T, nc)
    return log_activity_wj!(
        log_a,
        jacobian,
        aij,
        aux1,
        aux2,
        mix,
        nmol,
        volume,
        RT,
    )
end
