include("types.jl")
include("nvt.jl")
include("state_abstract.jl")
include("state_physical.jl")
include("state_idealidentity.jl")


"""
    vt_stability(mixture, nmol, volume, RT[, StateVariables; tol=1e-3, tpd_thresh=-1e-5, maxiter=200, buf])

Checks thermodynamical stability of single-phase state for `mixture` with composition `nmol`,
occupying `volume` at thermal energy `RT`.
Internally, runs direct minimization of dimensionless Helmholtz free energy TPD function (tangent plane distance)
from each of initial guess.

See also [`vt_stability!`](@ref).

*Optimizer.* BFGS with modified Cholesky decomposition (see `Downhill.CholBFGS`).

*Initialization.* Four guesses for the argument (Mikyska and Firoozabadi, Fluid Ph. Equil., 2012, 10.1016/j.fluid.2012.01.026).
Initial Hessian is put to the optimizer.

*`StateVariables.`* Available are
- [`VTStabilityPhysicalState`](@ref);
- [`VTStabilityIdealIdentityState`](@ref).

*Convergence criterion.*

```
  1
  -- maxᵢ (μ'ᵢ - μᵢ) < tolerance
  RT
```
where μ'ᵢ is chemical potential of trial phase and μᵢ is chemical potential of base (bulk) phase.

# Required arguments

- `mixture`: a mixture described by an equation of state;
- `nmol`: composition of the mixture, [mole];
- `volume::Real`: volume of the mixture [meter³];
- `RT`: thermal energy, `CubicEoS.GAS_CONSTANT_SI * temperature`, [Joule / mole];
- `StateVariables::Type{<:AbstractVTStabilityState}`: variables used in minimization,
    default is `VTStabilityPhysicalState`.

# Optional arguments

- `tol::Real=1e-3`: controls convergence criterion (see above);
- `tpd_thresh::Real=-1e-5`: negative TPD value, below which state is considered unstable, formally [moles];
- `maxiter::Int=200`: maximum number of iterations allowed per minimization;
- `buf::AbstractEoSThermoBuffer`: cache structure, reduces allocations, by default is `thermo_buffer(mixture)`.
"""
function vt_stability(mixture, nmol, volume, RT, ::Type{StateVariables};
    tol::Real=1e-3,
    tpd_thresh=-1e-5,
    maxiter::Int=1000,
    buf::AbstractEoSThermoBuffer=thermo_buffer(mixture),
) where {StateVariables<:AbstractVTStabilityState}
    # form base state
    basestate = VTStabilityBaseState(mixture, nmol, volume, RT; buf=buf)

    # prepare TPD closures: TPD, gradient, constrain_step
    tpd_fdf!, tpd_df! = __vt_stability_tpd_closures(StateVariables, basestate; buf=buf)
    constrain_step = __vt_stability_step_closure(StateVariables, basestate.mixture.components.b)

    # prepare stop criterion closure
    convcond = __vt_stability_convergence_closure(StateVariables, basestate, tol; buf=buf)

    # prepare initial guesses
    trial_concentrations = vt_stability_initials_satpressure(mixture, nmol, RT; buf=buf)

    # run optimizer for each guess
    optmethod = Downhill.CholBFGS(basestate.logconcentration)
    results = map(trial_concentrations) do concentration
        trialstate = fromconcentration(StateVariables, concentration)
        return vt_stability!(trialstate, basestate, optmethod;
            tpd_fdf! = tpd_fdf!,
            constrain_step=constrain_step,
            convcond=convcond,
            tpd_thresh=tpd_thresh,
            maxiter=maxiter,
            buf=buf,
        )
    end

    issuccess = any(x -> x.issuccess, results)
    !issuccess && error("VTStability: all tries have failed")

    # Not successed results are ignored
    isstable = all(x -> x.issuccess ? x.isstable : true, results)

    return issuccess, isstable, results
end

# Default variables
vt_stability(mixture, nmol, volume, RT; kwargs...) =
    vt_stability(mixture, nmol, volume, RT, VTStabilityPhysicalState; kwargs...)

function vt_stability!(
    trialstate::AbstractVTStabilityState,
    basestate::VTStabilityBaseState,
    optmethod;
    tpd_fdf!::Function,
    constrain_step::Function,
    convcond::Function,
    tpd_thresh::Real=-1e-5,
    maxiter::Int=200,
    buf::AbstractEoSThermoBuffer=thermo_buffer(basestate.mixture),
)
    testhessian = let n = ncomponents(basestate.mixture)
        # Matrix{Float64}(undef, (n, n))
        diagm(n, n, 0 => ones(n))
    end
    # If used hessian is true hessian
    # testhessian = helmholtztpdhessian!(testhessian, trialstate, basestate; buf=buf)

    teststatex = value(trialstate)
    Downhill.reset!(optmethod, teststatex, testhessian)

    isfinished = false
    optim = OptimStats(false, -1, -1)
    try
        optimresult = Downhill.optimize!(tpd_fdf!, optmethod, teststatex,
            convcond=convcond,
            maxiter=maxiter,
            constrain_step=constrain_step,
            reset=false,
        )
        isfinished = true
        optim = OptimStats(optimresult.converged, optimresult.iterations, optimresult.calls)
    catch e
        if e isa InterruptException
            rethrow(e)
        else
            @warn "$(sprint(showerror, e))"
        end
    end
    # Manually update trialstate
    teststatex .= Downhill.argumentvec(optmethod)

    # Result
    tpd_val = Downhill.fnval(optmethod)
    isstable = tpd_val ≥ -abs(tpd_thresh)
    testconc = concentration(trialstate)

    return VTStabilityResult(isfinished, isstable, tpd_val, testconc, trialstate, optim)
end

struct OptimStats
    converged::Bool
    iters::Int
    calls::Int
end

struct VTStabilityResult{T<:Real,S<:AbstractVTStabilityState}
    issuccess::Bool
    isstable::Bool
    energy_density::T
    concentration::Vector{T}
    state::S
    optim::OptimStats
end

function __vt_stability_tpd_closures(
    StateVariables::Type{<:AbstractVTStabilityState},
    basestate::VTStabilityBaseState;
    buf::AbstractEoSThermoBuffer=thermo_buffer(basestate.mixture),
)
    # for internal usage
    trialstate = StateVariables(similar(basestate.logconcentration))
    trialx = value(trialstate)

    function clsr_tpd_fdf!(x::AbstractVector, g::AbstractVector)
        trialx .= x
        # TODO: Avoid double calculation of gradient (`helmholtztpd` call)
        g = helmholtztpdgradient!(g, trialstate, basestate; buf=buf)
        tpd = helmholtztpd(trialstate, basestate; buf=buf)
        return tpd, g
    end

    function clsr_tpd_df!(g::AbstractVector, x::AbstractVector)
        trialx .= x
        g = helmholtztpdgradient!(g, trialstate, basestate; buf=buf)
        return g
    end

    return clsr_tpd_fdf!, clsr_tpd_df!
end

function __vt_stability_step_closure(
    ::Type{StateVariables},
    covolumes::AbstractVector,
) where {StateVariables<:AbstractVTStabilityState}
    function clsr(x::AbstractVector, direction::AbstractVector)
        return __constrain_step(StateVariables, x, direction, covolumes)
    end
    return clsr
end

function __vt_stability_convergence_closure(
    ::Type{StateVariables},
    basestate::VTStabilityBaseState,
    tol::Real;
    buf::AbstractEoSThermoBuffer=thermo_buffer(basestate.mixture),
) where {StateVariables<:AbstractVTStabilityState}
    trial_concentration = similar(basestate.logconcentration)
    trial_logcactivity = similar(basestate.logcactivity)
    chemdiff = similar(basestate.logcactivity)
    base_chempot = basestate.logconcentration + basestate.logcactivity

    function clsr_convcond(x::V, xpre::V, y::T, ypre::T, g::V) where {T<:Real,V<:AbstractVector}
        # No assignment here because of mutable capture
        concentration!(trial_concentration, StateVariables, x)
        log_c_activity!(
            trial_logcactivity,
            basestate.mixture,
            trial_concentration,
            1,
            basestate.RT;
            buf=buf,
        )
        @. chemdiff = log(trial_concentration) + trial_logcactivity - base_chempot
        return maximum(abs, chemdiff) < abs(tol)
    end
    return clsr_convcond
end

# TODO: Check if zfactors are distinct
# TODO: DRY
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
