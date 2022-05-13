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
