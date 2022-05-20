"EoS constraint with no bounds on magnitude of step for vt stability."
unconstrained_eos_step(::Type{<:AbstractVTStabilityState}, x, d) = (-Inf, Inf)

# TODO: Check if zfactors are distinct
# TODO: DRY
"Stability initial concentrations from ideal mixing rules (4 guesses)."
function vt_stability_initials_satpressure(
    mixture::AbstractEoSMixture,
    nmolbase::AbstractVector,
    RT::Real;
    buf::AbstractEoSThermoBuffer=thermo_buffer(mixture),
)
    psat = wilson_saturation_pressure.(mixture.components, RT)

    # Base is liquid, trial is vapor/liquid
    p_baseliquid = dot(psat, nmolbase) / sum(nmolbase)
    molfrac_trialvapor = (psat ./ p_baseliquid) .* (nmolbase ./ sum(nmolbase))
    zroots_baseliquid = zfactors(mixture, molfrac_trialvapor, p_baseliquid, RT; buf=buf)

    conc_trial_lv = let z = zfactorchoose(zroots_baseliquid, 'g')
        p_baseliquid .* molfrac_trialvapor ./ (z * RT)
    end
    conc_trial_ll = let z = zfactorchoose(zroots_baseliquid, 'l')
        p_baseliquid .* molfrac_trialvapor ./ (z * RT)
    end

    # Base is vapor, trial is vapor/liquid
    molfrac_trialliquid = (nmolbase ./ psat) ./ sum(nmolbase ./ psat)
    p_basevapor = dot(psat, molfrac_trialliquid)
    zroots_basevapor = zfactors(mixture, molfrac_trialliquid, p_basevapor, RT; buf=buf)

    conc_trial_vv = let z = zfactorchoose(zroots_basevapor, 'g')
        p_basevapor .* molfrac_trialliquid ./ (z * RT)
    end
    conc_trial_vl = let z = zfactorchoose(zroots_basevapor, 'l')
        p_basevapor .* molfrac_trialliquid ./ (z * RT)
    end

    return (conc_trial_lv, conc_trial_ll, conc_trial_vv, conc_trial_vl)
end
