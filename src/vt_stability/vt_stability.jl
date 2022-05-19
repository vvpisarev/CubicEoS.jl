include("types.jl")
include("nvt.jl")
include("state_abstract.jl")
include("state_physical.jl")
include("state_idealidentity.jl")
include("utils.jl")

"""
    vt_stability(mixture, nmol, volume, RT[, StateVariables; trial_concentrations, eos_constrain_step, tol=1e-3, tpd_thresh=-1e-5, maxiter=200, buf])

Checks thermodynamical stability of single-phase state for `mixture` with composition `nmol`,
occupying `volume` at thermal energy `RT`.
Internally, runs direct minimization of dimensionless Helmholtz free energy TPD function (tangent plane distance)
from each of initial guess.

See also [`vt_stability!`](@ref).

*Initialization.* Stability runs from a set of initial guesses provided by `trial_concentrations`.
By default, four guesses is used based on cubic eos and ideal mixing rules
(Mikyska and Firoozabadi, Fluid Ph. Equil., 2012, 10.1016/j.fluid.2012.01.026).
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

*Optimizer.* BFGS with modified Cholesky decomposition (see `Downhill.CholBFGS`).

*Optimization constraints.* There are two sets of constraints: physical and EoS-based.
Physical constraints are

```
  c'ᵢ > 0 for each component
```

and served by `CubicEoS` (see [`CubicEoS.physical_constrain_step`](@ref) for example).

EoS-based constraints are provided by user through `eos_constrain_step`.
For example, see [`CubicEoS.BrusilovskyEoS.eos_constrain_step`](@ref) which defines
"size of phase" constraints.
By default, there is no eos-related constraint.

# Required arguments

- `mixture`: a mixture described by an equation of state;
- `nmol`: composition of the mixture, [mole];
- `volume::Real`: volume of the mixture [meter³];
- `RT`: thermal energy, `CubicEoS.GAS_CONSTANT_SI * temperature`, [Joule / mole];
- `StateVariables::Type{<:AbstractVTStabilityState}`: variables used in minimization,
    default is `VTStabilityIdealIdentityState`.

# Optional arguments

- `trial_concentrations=nothing`: iterable of vectors of trial (initial) concentrations, by default four guesses from cubic eos is used (see 10.1016/j.fluid.2012.01.026);
- `eos_constrain_step::Function=(StateVariables, x, d) -> (-Inf, Inf)`: bounds for step magnitude `(αlow, αmax)` in optimization point `x` and direction `d` in variables `StateVariables` (`xnew = x + α*d`);
- `tol::Real=1e-3`: controls convergence criterion (see above);
- `tpd_thresh::Real=-1e-5`: negative TPD value, below which state is considered unstable, formally [moles];
- `maxiter::Int=1000`: maximum number of iterations allowed per minimization;
- `buf::AbstractEoSThermoBuffer`: cache structure, reduces allocations, by default is `thermo_buffer(mixture)`.

# Examples

Stability of methane (60%) + n-pentane (40%) at 300K and overall concentration 5000 mole / meter³

```julia
c1c5 = CubicEoS.load(BrusilovskyEoSMixture; names=("methane", "n-pentane"))
volume = 1e-6
RT = CubicEoS.GAS_CONSTANT_SI * 300
nmol = [0.6, 0.4] * 5000 * volume
vars = CubicEoS.VTStabilityIdealIdentityState

vt_stability(c1c5, nmol, volume, RT, vars;
    eos_constrain_step=CubicEoS.BrusilovskyEoS.eos_constrain_step(vars, c1c5),
)
```
"""
function vt_stability(mixture, nmol, volume, RT, ::Type{StateVariables};
    eos_constrain_step::Function=unconstrained_eos_step,
    trial_concentrations=nothing,
    tol::Real=1e-3,
    tpd_thresh=-1e-5,
    maxiter::Int=1000,
    buf::AbstractEoSThermoBuffer=thermo_buffer(mixture),
) where {StateVariables<:AbstractVTStabilityState}
    basestate = VTStabilityBaseState(mixture, nmol, volume, RT; buf=buf)

    tpd_fdf! = __vt_stability_tpd_closure(StateVariables, basestate)

    # Prepare stop criterion closure
    convcond = __vt_stability_convergence_closure(StateVariables, basestate, tol)

    # Prepare initial guesses, if not provided
    # Initital guess from a cubic eos
    if isnothing(trial_concentrations)
        trial_concentrations = let mixtureBrus = load(BrusilovskyEoSMixture; names=map(name, components(mixture)))
            vt_stability_initials_satpressure(mixtureBrus, nmol, RT)
        end
    end

    # Run optimizer for each guess
    optmethod = Downhill.CholBFGS(basestate.logconcentration)
    results = map(trial_concentrations) do concentration
        trialstate = fromconcentration(StateVariables, concentration)
        return vt_stability!(trialstate, basestate, optmethod;
            tpd_fdf! = tpd_fdf!,
            eos_constrain_step=eos_constrain_step,
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
    vt_stability(mixture, nmol, volume, RT, VTStabilityIdealIdentityState; kwargs...)

function vt_stability!(
    trialstate::StateVariables,
    basestate::VTStabilityBaseState,
    optmethod;
    tpd_fdf!::Function=__vt_stability_tpd_closure(StateVariables, basestate),
    eos_constrain_step::Function=unconstrained_eos_step,
    convcond::Function,
    tpd_thresh::Real=-1e-5,
    maxiter::Int=200,
    buf::AbstractEoSThermoBuffer=thermo_buffer(basestate.mixture),
) where {StateVariables<:AbstractVTStabilityState}
    constrain_step = __vt_stability_step_closure(StateVariables, eos_constrain_step)

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

function __vt_stability_tpd_closure(
    StateVariables::Type{<:AbstractVTStabilityState},
    basestate::VTStabilityBaseState;
    buf::AbstractEoSThermoBuffer=thermo_buffer(basestate.mixture),
)
    trialstate = StateVariables(similar(basestate.logconcentration))
    trialx = value(trialstate)

    function clsr_tpd_fdf!(x::AbstractVector, g::AbstractVector)
        trialx .= x
        # TODO: Avoid double calculation of gradient (`helmholtztpd` call)
        g = helmholtztpdgradient!(g, trialstate, basestate; buf=buf)
        tpd = helmholtztpd(trialstate, basestate; buf=buf)
        return tpd, g
    end

    return clsr_tpd_fdf!
end

"Creates closure-constraint for optimization."
function __vt_stability_step_closure(
    ::Type{StateVariables},
    eos_constrain_step::Function,
) where {StateVariables<:AbstractVTStabilityState}
    function clsr(x::AbstractVector, direction::AbstractVector)
        αlowphys, αmaxphys = physical_constrain_step(StateVariables, x, direction)
        αloweos, αmaxeos = eos_constrain_step(StateVariables, x, direction)

        αlow = max(αlowphys, αloweos)
        αmax = min(αmaxphys, αmaxeos)

        (αmax ≤ 0 || αmax ≤ αlow) && throw(ConstrainStepLowerBoundError(x, direction))

        return αmax
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
