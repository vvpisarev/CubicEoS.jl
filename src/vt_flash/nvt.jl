"""
    nvtgradient!(grad, mix, nmol, volume, RT[; buf])

Isothermal `grad`ient of dimensionless [Joules / RT] Helmholtz free energy of `mix`ture
at `nmol`, `volume` and `RT`.

`grad`ient should be column of size nc × 1, where nc is number of components.
To avoid intermediate allocations, use `buf`, see [`thermo_buffer`](@ref).

# Gradient's structure

```
1
-- [ μ₁, ..., μₙ, -P ]
RT
```
where μ is chemical potential and P is pressure.

See also [`nvtgradienthessian!`](@ref).
"""
function nvtgradient!(
    grad::AbstractVector,
    mix::BrusilovskyEoSMixture,
    nmol::AbstractVector,
    volume::Real,
    RT::Real;
    buf::BrusilovskyThermoBuffer=thermo_buffer(mix),
)
    # Partial derivatives by moles, ln Nᵢ - ln V + ln cₐᵢ
    dN = @view grad[1:end-1]
    dN = log_c_activity!(dN, mix, nmol, volume, RT; buf=buf)
    dN .+= log.(nmol)
    dN .-= log(volume)

    # Partial derivative by volume, -P/RT
    grad[end] = - pressure(mix, nmol, volume, RT; buf=buf) / RT

    return grad
end

"""
    nvthessian!(hess, mix, nmol, volume, RT[; buf])

Isothermal Hessian of dimensionless [Joules / RT] Helmholtz free energy of `mix`ture
at `nmol`, `volume` and `RT`.

`hess`ian should be matrix of size (nc + 1) × (nc + 1).
To avoid intermediate allocations, use `buf`, see [`thermo_buffer`](@ref).

# Hessian's structure

```
| B | C |  B is n×n, ∂²a/∂Nᵢ∂Nⱼ
| --+-- |  C is column n×1, ∂²a/∂Nᵢ∂V
| Cᵀ| D |  D is number 1×1, ∂²a/∂V²
```

See also [`nvtgradienthessian!`](@ref).
"""
function nvthessian!(
    hess::AbstractMatrix,
    mix::BrusilovskyEoSMixture,
    nmol::AbstractVector,
    volume::Real,
    RT::Real;
    buf::BrusilovskyThermoBuffer=thermo_buffer(mix),
)
    # Filling B block. Part [C; D] is used as auxiliary
    aux = @view hess[1:end-1, end]
    dNdN = @view hess[1:end-1, 1:end-1]
    aux, dNdN = log_c_activity_wj!(aux, dNdN, mix, nmol, volume, RT; buf=buf)
    # Adding diagonal term
    @inbounds for i in eachindex(nmol)
        dNdN[i, i] += 1 / nmol[i]
    end

    # Filling [C; D] blocks. It is gradient of pressure, actually
    dP = @view hess[1:end, end]
    # dP = [∂P/∂Nᵢ... ∂P/∂V]
    dP = vtpressuregradient!(dP, mix, nmol, volume, RT; buf=buf)
    dP ./= -RT  # scale by -1/RT
    # Copy C to Cᵀ
    hess[end, 1:end-1] .= hess[1:end-1, end]

    return hess
end

"""
    nvtgradienthessian!(grad, hess, mix, nmol, volume, RT[; buf])

Isothermal gradient and hessian of dimensionless [Joules / RT] Helmholtz free energy
of `mix`ture at `nmol`, `volume` and `RT`.

This version should be faster than
pair of calls [`nvtgradient!`](@ref), [`nvthessian!`](@ref).

For `grad` and `hess` see [`nvtgradient!`](@ref), [`nvthessian!`](@ref).
To avoid intermediate allocations, use `buf`, see [`thermo_buffer`](@ref).
"""
function nvtgradienthessian!(
    grad::AbstractVector,
    hess::AbstractMatrix,
    mix::BrusilovskyEoSMixture,
    nmol::AbstractVector,
    volume::Real,
    RT::Real;
    buf::BrusilovskyThermoBuffer=thermo_buffer(mix),
)
    # Hessian
    # `grad` is used here as auxiliary

    ## Second partial derivatives by (nmolᵢ and volume)
    ## Simultaneously, second partial derivate by volume
    dP = @view hess[1:end, end]
    dP = vtpressuregradient!(dP, mix, nmol, volume, RT; buf=buf)

    hess[1:end, end] ./= -RT
    hess[end, 1:end-1] .= hess[1:end-1, end]  # Symmetric

    ## Second partial derivatives by (nmolᵢ and nmolⱼ)
    dNdN = @view hess[1:end-1, 1:end-1]
    logca = @view grad[1:end-1]
    logca, dNdN = log_c_activity_wj!(logca, dNdN, mix, nmol, volume, RT; buf=buf)
    ## Adding diagonal term
    @inbounds for i in eachindex(nmol)
        dNdN[i, i] += 1 / nmol[i]
    end

    # Gradient

    ## Partial derivatives by moles, ln Nᵢ - ln V + ln cₐᵢ
    dN = logca
    dN .+= log.(nmol)
    dN .-= log(volume)

    ## Partial derivative by volume, -P/RT
    grad[end] = -pressure(mix, nmol, volume, RT; buf=buf) / RT

    return grad, hess
end

"""
    vtpressuregradient!(grad, mix, nmol, volume, RT[, buf])

Pressure `grad`ient at point `[nmol...; volume]` and `RT`.
To avoid intermediate allocations, use `buf`, see [`thermo_buffer`](@ref).

*The gradient does not include ∂P/∂T.*

# Gradient's structure

```
[∂P/∂N₁, ..., ∂P/∂Nₙ, ∂P/∂V]
```

where `n` is number of components in `mix`ture.
"""
function vtpressuregradient!(
    ∇P::AbstractVector{T},
    mix::BrusilovskyEoSMixture{T},
    nmol::AbstractVector,
    volume::Real,
    RT::Real;
    buf::BrusilovskyThermoBuffer=thermo_buffer(mix),
) where {T}
    # I did not implement this function in src/basic_thermo.jl
    # because the gradient here does not include ∂P/∂T derivative.
    # Maybe, it should be implemented later in src/basic_thermo.jl.

    A, B, C, D, aij = eos_parameters(mix, nmol, RT; buf=buf)

    # hell arithmetics
    # does compiler smart to detect this as constants
    # if plain operation were put in ∂P/∂Nᵢ for-cycle explicitly?
    V = volume  # alias
    VmB⁻¹ = 1 / (V - B)
    ΣnmolbyVmB² = sum(nmol) * VmB⁻¹^2
    DmC = D - C
    VpC⁻¹ = 1 / (V + C)
    VpC⁻² = VpC⁻¹^2
    VpD⁻¹ = 1 / (V + D)
    VpD⁻² = VpD⁻¹^2
    AbyDmC = A / DmC
    VpC⁻¹mVpD⁻¹byDmC² = (VpC⁻¹ - VpD⁻¹) / DmC^2

    # ∂P/∂Nᵢ part
    for (i, substance) in enumerate(components(mix))
        bᵢ, cᵢ, dᵢ = substance.b, substance.c, substance.d
        ∂ᵢA = 2 * dot(nmol, @view aij[i, :])  # ∂A/∂Nᵢ

        ∇P[i] = RT * (VmB⁻¹ + bᵢ * ΣnmolbyVmB²) - (
            (∂ᵢA * DmC - A * (dᵢ - cᵢ)) * VpC⁻¹mVpD⁻¹byDmC²
            + AbyDmC * (-cᵢ * VpC⁻² + dᵢ * VpD⁻²)
        )
    end
    ∇P[end] = - RT * ΣnmolbyVmB² + AbyDmC * (VpC⁻² - VpD⁻²)
    return ∇P
end
