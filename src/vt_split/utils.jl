"EoS constraint with no bounds on magnitude of step for vt phase split."
unconstrained_eos_step(::Type{<:AbstractVTFlashState}, x, d) = (-Inf, Inf)

"""
    __sort_phases!(mix, nmol₁, V₁, nmol₂, V₂, RT) -> (nmolgas, volumegas, nmolliq, volumeliq)

Sort moles and volumes of phases 1 and 2 into gas and liquid by their z-factors.
Values in `nmol₁` and `nmol₂` vectors may be swapped.
"""
function __sort_phases!(mix, nmol₁, V₁, nmol₂, V₂, RT)
    P₁ = pressure(mix, nmol₁, V₁, RT)  # Pressures should be equal
    P₂ = pressure(mix, nmol₂, V₂, RT)

    # TODO: Seems can be reduced to Vᵢ / sum(nmolᵢ)
    Z₁ = P₁ * V₁ / (sum(nmol₁) * RT)
    Z₂ = P₂ * V₂ / (sum(nmol₂) * RT)

    if Z₂ > Z₁  # □₂ is gas state, need exchange
        P₁, P₂ = P₂, P₁
        Z₁, Z₂ = Z₂, Z₁
        V₁, V₂ = V₂, V₁

        for i in eachindex(nmol₁, nmol₂)
            nmol₁[i], nmol₂[i] = nmol₂[i], nmol₁[i]
        end
    end
    # now □₁ is gas, □₂ is liquid
    return nmol₁, V₁, nmol₂, V₂
end

# TODO: indexes may be inconsistent with enumerate
# TODO: may be implement with `minimum(f, A)` (and `filter`) function?
"""
Return concentration of state with minimum energy from iterable of vt-stability tries.
Formely used in flash as default strategy for choice of initial concentration.
"""
function concentrationwithlowesttpd(vt_stab_tries)
    Dmin = Inf
    index_min = -1
    for (i, state) in enumerate(vt_stab_tries)
        if state.issuccess && !state.isstable && state.energy_density < Dmin
            index_min = i
            Dmin = state.energy_density
        end
    end
    index_min == -1 && error("Stability tries are inconsistent: Can't choose the one with the lowest energy")
    return vt_stab_tries[index_min].concentration
end

vtsplitvariables() = (VTFlashPhysicalState, VTFlashRatioState, VTFlashIdealIdentityState)
