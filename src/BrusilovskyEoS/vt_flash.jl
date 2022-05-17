function CubicEoS.vtpressuregradient!(
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

# General-purpose covolume constraint based on bisection and moles-volume variables.
function eos_vt_split_constrain_step(
    StateVariables::Type{<:CubicEoS.AbstractVTFlashState},
    mix::BrusilovskyEoSMixture{T},
    nmolbase::AbstractVector,
    volumebase::Real,
) where {T}
    maxiter = 50

    xtrial = Vector{T}(undef, ncomponents(mix) + 1)
    covolumes = mix.components.b
    nmol2 = similar(nmolbase, T)
    nmoltrial = similar(nmolbase, T)

    function iscovmeet(nmol1, volume1)
        @. nmol2 = nmolbase - nmol1
        volume2 = volumebase - volume1

        cov1 = dot(nmol1, covolumes) ≤ volume1
        cov2 = dot(nmol2, covolumes) ≤ volume2
        return cov1 && cov2
    end

    function clsr(
        StateVariables::Type{<:CubicEoS.AbstractVTFlashState},
        x,
        direction,
    )
        αlow = zero(T)
        αmax = 1e3

        for _ in 1:maxiter
            @. xtrial = x + αmax * direction
            # If `_` is replaced with `nmoltrial`, type inference crashes
            _, volumetrial = CubicEoS.nmolvol!(nmoltrial, StateVariables(xtrial), nmolbase, volumebase)
            iscovmeet(nmoltrial, volumetrial) && break
            αmax *= 0.5
        end
        return αlow, αmax
    end
    return clsr
end

# Effective implementation of covolume constraint for physical variables.
function eos_vt_split_constrain_step(
    StateVariables::Type{<:CubicEoS.VTFlashPhysicalState},
    mix::BrusilovskyEoSMixture,
    nmolbase::AbstractVector,
    volumebase::Real,
)
    covolumes = [mix.components.b; -1.0]
    basedotcov = dot([nmolbase; volumebase], covolumes)

    function clsr(
        StateVariables::Type{<:CubicEoS.VTFlashPhysicalState},
        x,
        direction,
    )
        xdotcov = dot(x, covolumes)
        dirdotcov = dot(direction, covolumes)

        αmax = dirdotcov > 0 ? - xdotcov / dirdotcov : (basedotcov - xdotcov) / dirdotcov
        αlow = dirdotcov < 0 ? - xdotcov / dirdotcov : (basedotcov - xdotcov) / dirdotcov

        return αlow, αmax
    end
end

# Effective implementation of covolume constraint for ratio variables.
# TODO: DRY: Seems very similar to implememntation for physical variables,
# TODO: DRY: only captured variables differ.
function eos_vt_split_constrain_step(
    StateVariables::Type{<:CubicEoS.VTFlashRatioState},
    mix::BrusilovskyEoSMixture,
    nmolbase::AbstractVector,
    volumebase::Real,
)
    covolumes = [nmolbase .* mix.components.b; -volumebase]
    basedotcov = sum(covolumes)

    function clsr(
        StateVariables::Type{<:CubicEoS.VTFlashRatioState},
        x,
        direction,
    )
        xdotcov = dot(x, covolumes)
        dirdotcov = dot(direction, covolumes)

        αmax = dirdotcov > 0 ? - xdotcov / dirdotcov : (basedotcov - xdotcov) / dirdotcov
        αlow = dirdotcov < 0 ? - xdotcov / dirdotcov : (basedotcov - xdotcov) / dirdotcov

        return αlow, αmax
    end
end
