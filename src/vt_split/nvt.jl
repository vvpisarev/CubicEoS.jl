"""
    helmholtz(mix, nmol, volume, RT[; buf])

Dimensionless Helmholtz free energy [Joules / RT] of `mix`ture
at `nmol`, `volume` and `RT`.

To reduce intermediate allocations, use `buf`, see [`thermo_buffer`](@ref).
"""
function helmholtz(
    mix::AbstractEoSMixture,
    nmol::AbstractVector,
    volume::Real,
    RT::Real;
    buf::AbstractEoSThermoBuffer=thermo_buffer(mix),
)
    # helmholtz = dot([nmol; V]ᵀ, grad)

    grad = similar(nmol, Float64, length(nmol) + 1)
    grad = nvtgradient!(grad, mix, nmol, volume, RT; buf=buf)
    grad[1:end-1] .*= nmol
    grad[end] *= volume
    return sum(grad)
end

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
    mix::AbstractEoSMixture,
    nmol::AbstractVector,
    volume::Real,
    RT::Real;
    buf::AbstractEoSThermoBuffer=thermo_buffer(mix),
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
    mix::AbstractEoSMixture,
    nmol::AbstractVector,
    volume::Real,
    RT::Real;
    buf::AbstractEoSThermoBuffer=thermo_buffer(mix),
)
    # Filling B block. Part C is used as auxiliary
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
    mix::AbstractEoSMixture,
    nmol::AbstractVector,
    volume::Real,
    RT::Real;
    buf::AbstractEoSThermoBuffer=thermo_buffer(mix),
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

# Gradient's structure

```
[∂P/∂N₁, ..., ∂P/∂Nₙ, ∂P/∂V]
```

where `n` is number of components in `mix`ture.

*Note that, the gradient does not include ∂P/∂T.*
"""
function vtpressuregradient!(
    ∇P::AbstractVector{T},
    mix::AbstractEoSMixture{T},
    nmol::AbstractVector,
    volume::Real,
    RT::Real;
    buf::AbstractEoSThermoBuffer=thermo_buffer(mix),
) where {T}
    error("NotImplemented")
end
