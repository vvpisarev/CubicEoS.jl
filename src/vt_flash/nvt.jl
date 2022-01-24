"""
    nvtgradient!(grad, mix, nmol, volume, RT[; buf])

Gradient of dimensionless [Joules / RT] Helmholtz free energy at `nmol`, `volume` and `RT`.

To avoid intermediate allocations, use `buf`, see [`thermo_buffer`](@ref).

# Gradient's structure

```
1
-- [ μ₁, ..., μₙ, -P ]
RT
```
where μ is chemical potential and P is pressure.
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
    nvtgradienthessian!(grad, hess, mix, nmol, volume, RT[; buf])

Gradient and hessian of dimensionless [Joules / RT] Helmholtz free energy
at `nmol`, `volume` and `RT`.

To avoid intermediate allocations, use `buf`, see [`thermo_buffer`](@ref).

# Gradient's structure

See [`nvtgradient!`](@ref).

# Hessian's structure

```
| B | C |  B is n×n, ∂²a/∂Nᵢ∂Nⱼ
| --+-- |  C is column n×1, ∂²a/∂Nᵢ∂V
| Cᵀ| D |  D is number 1×1, ∂²a/∂V²
```
"""
function nvtgradienthessian!(
    grad::AbstractVector,  # n+1
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
    aux = grad
    aux = __vt_flash_pressure_gradient!(aux, mix, nmol, volume, RT; buf=buf)

    dPdN = @view aux[1:end-1]
    dPdV = aux[end]

    hess[1:end-1, end] .= .- dPdN ./ RT
    hess[end, 1:end-1] .= hess[1:end-1, end]  # Symmetric

    hess[end, end] = -1 * dPdV / RT  # d²a / dV²

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
    grad[end] = .- pressure(mix, nmol, volume, RT; buf=buf) / RT

    return grad, hess
end
