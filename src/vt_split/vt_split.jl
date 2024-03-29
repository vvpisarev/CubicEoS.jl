include("types.jl")
include("nvt.jl")
include("state_abstract.jl")
include("state_physical.jl")
include("state_ratio.jl")
include("state_idealidentity.jl")
# include("newton.jl")
include("utils.jl")

"""
    vt_split(mixture, nmol, volume, RT, trial_concetration, StateVariables[;
        tol, chemtol=tol, presstol=tol, maxiter=100, eos_constrain_step=unconstrained_eos_step]
    )

Two-phase thermodynamical phase split solver for `mixture` at given moles `nmol`, `volume`
and thermal energy `RT` from `trial_concentration` of a phase.
Phase split is performed with `StateVariables` using Cholesky-BFGS optimization.

The phase equilibrium is considered, when

1. Chemical potentials are equal in a sense

```
 1
--- maxᵢ |μᵢ' - μᵢ''| < chemtol
 RT
```

2. Pressures are equals in a sense
```
|P' - P''| volume
----------------- < presstol,
   RT sum(nmol)
```
where `i` is component index, and `'`, `''` are indexes of phases.

Return [`VTSplitResult`](@ref).

See also [`CubicEoS.vt_split!`](@ref), [`vt_flash_newton`](@ref).

# Arguments

- `mix::AbstractEoSMixture{T}`: mixture;
- `nmol::AbstractVector`: moles of mixtures' components [mole];
- `volume::Real`: volume of mixture [meter³];
- `RT::Real`: thermal energy in [Joule / mole], to specify from temperature use `CubicEoS.GAS_CONSTANT_SI * temperature`;
- `trial_concentration::AbstractVector`: concentration of components in a trial state, should corresponds to negative Helmholtz TPD;
- `StateVariables::Type{<:AbstractVTSplitState}`: one of state variables, call `CubicEoS.vtsplitvariables()` to list.

# Optional arguments

- `chemtol::Real=tol`: tolerance for chemical potentials of components;
- `presstol::Real=tol`: tolerance for pressures of phases;
- `tol::Real=1024*eps(T)`: default tolerance for both `chemtol` and `presstol`,
    `T` is AbstractFloat type defined by `mixture`'s type;
- `maxiter::Integer`: maximum allowed steps in optimization;
- `eos_constrain_step::Function=CubicEoS.unconstrained_eos_step`:
    function of signature `(StateVariables, x, d) -> (αlow, αmax)`.
    If you have a constraint from EoS of form `A ≤ x + α d ≤ B`,
    this function returns bounds `αlow` and `αmax` on magnitude of step `α` in optimization step,
    `x` is `StateVariables` value in the step,
    `d` is direction (in same variables).
"""
function vt_split(
    mix::AbstractEoSMixture{T},
    nmol::AbstractVector,
    volume::Real,
    RT::Real,
    concentration::AbstractVector,
    StateVariables::Type{<:AbstractVTSplitState};
    tol::Real=1024*eps(T),
    chemtol::Real=tol,
    presstol::Real=tol,
    maxiter::Integer=100,
    eos_constrain_step::Function=unconstrained_eos_step,
) where {T}
    # Renaming
    saturation = __find_saturation_negative_helmdiff(mix, nmol, volume, RT, concentration;
        maxsaturation=1,
        maxiter=50,
        scale=0.5,
        helmdifftol=-1e-7/RT,
    )

    if isnan(saturation)
        @error "VTSplit: Initial state was not found!" mixture=mix nmol=repr(nmol) volume=volume RT=RT
        error("VTSplit: Initial state was not found!")
    end

    state = StateVariables(concentration, saturation, nmol, volume)

    return vt_split!(
        state,
        mix,
        nmol,
        volume,
        RT;
        chemtol=chemtol,
        presstol=presstol,
        maxiter=maxiter,
        eos_constrain_step=eos_constrain_step,
    )
end

"""
    vt_split!(unstable_state, mix, nmol, volume, RT; chemtol, presstol, maxiter)

VT phase split from an `unstable_state::AbstractVTSplitState`, which will be destructed.

For rest of arguments see [`vt_split`](@ref).

Return [`VTSplitResult`](@ref).
"""
function vt_split!(
    unstable_state::StateVariables,
    mix::AbstractEoSMixture,
    nmol::AbstractVector,
    volume::Real,
    RT::Real;
    chemtol::Real,
    presstol::Real,
    maxiter::Int,
    eos_constrain_step::Function=unconstrained_eos_step,
) where {StateVariables<:AbstractVTSplitState}
    state = unstable_state
    statex = value(state)

    # Initial hessian
    hessian = Matrix{Float64}(undef, (size(statex, 1), size(statex, 1)))
    hessian = hessian!(hessian, state, mix, nmol, volume, RT)

    # hessian = diagm(0 => ones(size(statex, 1)))

    # Construction of step-closure from physical and eos constraints
    constrain_step = let
        physical_constrain_step = physical_constrain_step_closure(
            StateVariables,
            physical_constrain_step_lowlims(StateVariables, nmol, volume),
            physical_constrain_step_uplims(StateVariables, nmol, volume),
        )
        __vt_split_step_closure(
            StateVariables,
            physical_constrain_step,
            eos_constrain_step,
        )
    end

    helmdiff! = __vt_split_helmdiff_closure(state, mix, nmol, volume, RT)

    convcond = __vt_split_convergence_closure(state, mix, nmol, volume, RT;
        chemtol=chemtol,
        presstol=presstol,
    )

    # Run optimizer
    optmethod = Downhill.CholBFGS(statex)
    Downhill.reset!(optmethod, statex, hessian)
    optimresult = Downhill.optimize!(helmdiff!, optmethod, statex;
        gtol=NaN,
        convcond=convcond,
        maxiter=maxiter,
        constrain_step=constrain_step,
        reset=false,
    )
    statex .= optimresult.argument

    return VTSplitResult(state, mix, nmol, volume, RT, optimresult)
end

"Creates closure of helmholtz energy for optimization."
function __vt_split_helmdiff_closure(
    state1::AbstractVTSplitState,
    mix::AbstractEoSMixture,
    nmolb::AbstractVector,  # base state
    volumeb::Real,  # base state
    RT::Real,
)
    state1x = value(state1)  # link to internal
    buf = thermo_buffer(mix)

    helmb = helmholtz(mix, nmolb, volumeb, RT; buf=buf)

    nmol2 = similar(nmolb, Float64)

    "Takes `x`, modifies gradient, `state` is kept fresh."
    function clsr!(x::AbstractVector, grad::AbstractVector)
        state1x .= x

        nmol1, vol1 = nmolvol(state1, nmolb, volumeb)
        a1 = helmholtz(mix, nmol1, vol1, RT; buf=buf)

        @. nmol2 = nmolb - nmol1
        vol2 = volumeb - vol1
        a2 = helmholtz(mix, nmol2, vol2, RT; buf=buf)

        Δa = (a1 + a2) - helmb
        grad = gradient!(grad, state1, mix, nmolb, volumeb, RT; buf=buf)
        return Δa, grad
    end

    return clsr!
end

# TODO: DRY: same code as for `__vt_stability_step_closure`
function __vt_split_step_closure(
    ::Type{StateVariables},
    physical_constrain_step::Function,
    eos_constrain_step::Function,
) where {StateVariables<:AbstractVTSplitState}
    function clsr(x::AbstractVector, direction::AbstractVector)
        # Determine bounds from physical and eos constraints
        # Merge the bounds
        # Check for consistency of maximum step
        αlowphys, αmaxphys = physical_constrain_step(StateVariables, x, direction)
        αloweos, αmaxeos = eos_constrain_step(StateVariables, x, direction)

        αlow = max(αlowphys, αloweos)
        αmax = min(αmaxphys, αmaxeos)

        (αmax ≤ 0 || αmax ≤ αlow) && throw(ConstrainStepLowerBoundError(x, direction))

        return αmax
    end
    return clsr
end

# TODO: check that saturation is valid for physical and eos constraints
function __find_saturation_negative_helmdiff(
    mix::AbstractEoSMixture,
    nmolb::AbstractVector,
    volumeb::Real,
    RT::Real,
    concentration::AbstractVector;
    maxsaturation::Real,
    maxiter::Int,
    scale::Real,
    helmdifftol::Real,
    thermo_buf::AbstractEoSThermoBuffer=thermo_buffer(mix),
)
    # Δa = a' + a'' - abase
    abase = helmholtz(mix, nmolb, volumeb, RT; buf=thermo_buf)

    nmol1 = similar(nmolb, Float64)
    nmol2 = similar(nmolb, Float64)

    function helmdiff(saturation::Real)
        volume1 = volumeb * saturation
        @. nmol1 = volume1 * concentration
        a1 = helmholtz(mix, nmol1, volume1, RT; buf=thermo_buf)

        volume2 = volumeb - volume1
        @. nmol2 = nmolb - nmol1
        a2 = helmholtz(mix, nmol2, volume2, RT; buf=thermo_buf)
        return (a1 + a2) - abase
    end

    saturation = float(maxsaturation)
    for i in 1:maxiter
        try
            Δa = helmdiff(saturation)
            Δa < - abs(helmdifftol) && return saturation
        catch e
            isa(e, UndefVarError) && throw(e)  # Syntax
            # Other exceptions include invalid saturations, so, ignored
        end
        saturation *= scale
    end
    return NaN
end

function __vt_split_convergence_closure(
    state1::AbstractVTSplitState,
    mix::AbstractEoSMixture{T},
    nmolb::AbstractVector,
    volumeb::Real,
    RT::Real;
    chemtol::Real,
    presstol::Real,
    buf::AbstractEoSThermoBuffer=thermo_buffer(mix),
) where {T}
    state1x = value(state1)

    g1 = Vector{T}(undef, ncomponents(mix) + 1)
    g2 = similar(g1)
    diff = similar(g1)
    nmol2 = similar(g1, T, ncomponents(mix))

    function convcond(x::V, xpre::V, y::T1, ypre::T1, g::V) where {T1<:Real,V<:AbstractVector}
        state1x .= x
        nmol1, vol1 = nmolvol(state1, nmolb, volumeb)
        g1 = nvtgradient!(g1, mix, nmol1, vol1, RT; buf=buf)

        @. nmol2 = nmolb - nmol1
        vol2 = volumeb - vol1
        g2 = nvtgradient!(g2, mix, nmol2, vol2, RT; buf=buf)

        #=
                   [μ₁' - μ₁'', ..., μₙ' - μₙ''; -P' + P'']ᵀ
            diff = ----------------------------------------
                                    RT
        =#
        @. diff = g1 - g2

        condchem = let diffchem = (@view diff[1:end-1])
            # max_i |chem'_i - chem''_i| / RT < tolerance
            maximum(abs, diffchem) < chemtol
        end
        condpress = let diffpress = diff[end]
            # |P' - P''| V
            # ------------ < tolerance
            #    RT ΣNᵢ
            abs(diffpress) * volumeb / sum(nmolb) < presstol
        end

        return condchem && condpress
    end

    return convcond
end

"""
Extracts vt-state from `optresult` (Downhill obj).
Sorts variables into gas and liquid.
Returns corresponding `VTSplitResult`.
"""
function VTSplitResult(
    state::S,
    mix::AbstractEoSMixture{T},
    nmol::AbstractVector{T},
    volume::Real,
    RT::Real,
    optresult,
) where {T, S<:AbstractVTSplitState}
    nmol1, volume1 = nmolvol(state, nmol, volume)

    nmol2 = nmol .- nmol1
    volume2 = volume - volume1

    nmolgas, volgas, nmolliq, volliq = __sort_phases!(mix, nmol1, volume1, nmol2, volume2, RT)

    return VTSplitResult{T, S}(;
            converged=optresult.converged,
            singlephase=false,
            RT=RT,
            nmolgas=nmolgas,
            volumegas=volgas,
            nmolliq=nmolliq,
            volumeliq=volliq,
            state=state,
            iters=optresult.iterations,
            calls=optresult.calls,
    )
end
