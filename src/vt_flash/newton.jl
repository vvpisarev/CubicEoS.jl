using LinearAlgebra

"""
Result of Newton's iterations.
Contain `converged::Bool` flag, the minimum `argument::Vector`,
number of `iterations::Integer` and `calls::Integer`.
"""
struct NewtonResult{T}
    converged::Bool
    argument::Vector{T}
    iterations::Int
    calls::Int

    NewtonResult(conv, x, iters=-1, calls=-1) = new{Float64}(conv, copy(x), iters, calls)
end

"""
    backtracking_line_search(f, x₀, d, y₀[; α₀, p, buf])

Find `x = x₀ + α*d` (`x₀::AbstractVector`, `d::AbstractVector`) that `f(x)::Real`
is < `y₀ == f(x₀)` by exponentially decreasing `::Real` `α = α₀ * pⁱ` (`p::Real=0.5` must be < 1).
Optional `buf` is similar to `x₀` for storing `x₀ + α*d` trial state.

Return tuple of `α`, `f(x₀ + α*d)` and number of `f`'s calls.
"""
function backtracking_line_search(
    f::Function,
    x₀::AbstractVector{T},
    d::AbstractVector{T},
    y₀::T;
    α::T=one(T),
    p::Real=one(T)/2,
    buf::AbstractVector{T}=similar(x₀),
) where {T<:Real}
    xtry = buf
    calls = 0
    ftry = T(NaN)

    # (?) better `maxiter = 1 + ceil(Int, log(eps(α))/log(p))`
    # If use above, then `p^maxiter ≈ eps(α₀)`
    maxiter = 200
    for i in 1:maxiter
        @. xtry = x₀ + α*d
        calls += 1
        ftry = f(xtry)

        @debug "backtracking_line_search" i α p ftry y₀ ftry-y₀

        ftry < y₀ && break
        α *= p

        i == maxiter && error("Line search number of iterations ($(maxiter)) exceeded.")
    end
    return α, ftry, calls
end

"""
    newton(f, grad!, hess!, x₀[; gtol=1e-6, convcond, maxiter=200, constrain_step])

Find (un)constrained minimum of `f(x)::Real` using Newton's solver with line search
and hessian modified by Cholesky factorization.
The optimization starts at `x₀` performing no more than `maxiter` Newton's steps.

New point `x ← xpre + αmax * d` is constrained in Newton's `d`irection by
`constrain_step(x, d) -> α::Real`.
By default, the constrains are disabled.

Return [`NewtonResult`](@ref) object.

# Arguments

The target function and initialization of optimization is defined by following arguments

- `f(x) -> Real`: function to optimize;
- `grad!(g, x) -> g`: must update and return `g`radient of `f` at `x`;
- `hess!(H, x) -> H`: must update and return `H`essian matrix of `f` at `x`.
- `x₀::AbstractVector`: initial argument for optimization.

# Optional arguments

## Convergence criteria

By default, the optimization considered converged when `norm(∇(x), 2) ≤ gtol`.
Generally, the convergence criterion is defined by `convcond`.

- `gtol::Real=1e-6`: (default criterion) 2-norm of gradient when iterations
    are considered converged;
- `convcond::Function` with signature `(x, xpre, y, ypre, g)->Bool`:
    if returns `true` for given arguments `x` and `xpre`, function values `y=f(x)`,
    `ypre=f(xpre)` and gradient `g=∇(x)`, the optimization is stopped.

Example (default criterion): `stopcond=(x, xpre, y, ypre, g)->norm(g, 2) ≤ gtol`.

## Defining constrains

- `constrain_step::Function=(x, d)->Inf`:
    given an argument `x` and `d`irection, returns maximum allowed scalar `α` in form of
    `xnew = x + α*d`. Default is unconstrained optimization (`α = Inf`).

## Forcing stop of optimization

- `maxiter::Integer=200`: maximum allowed Newton's steps.
"""
function newton(
    f::Function,
    ∇f!::Function,
    H!::Function,
    x_::AbstractVector;
    gtol::Real=1e-6,
    maxiter::Integer=200,
    constrain_step::Function=(x, d)->Inf,
    convcond::Function=(x, xpre, y, ypre, g)->norm(g, 2) ≤ gtol,
)
    x = float.(x_)  # copy x for internal use

    ∇ = similar(x)
    δx = similar(x)
    hess_full = Matrix{eltype(x)}(undef, size(∇, 1), size(∇, 1))
    vec = similar(x)
    xpre = similar(x)
    totfcalls = 0


    """
    Returns closure for Downhill.strong_backtracking!
    Vector `v` for determine types.
    """
    function fdfclosure(v::AbstractVector)
        vecx = similar(v)
        vecg = similar(v)

        function fdf(x::AbstractVector, α::Real, d::AbstractVector)
            vecx .= x .+ α .* d
            return f(vecx), ∇f!(vecg, vecx)
        end

        return fdf
    end

    fdf = fdfclosure(x)

    # Initializing, main cycle is postconditioned
    y = f(x)
    totfcalls += 1
    ∇ = ∇f!(∇, x)

    for i in 1:maxiter
        # descent direction `δx`

        hess_full = H!(hess_full, x)
        hess = Downhill.mcholesky!(hess_full)

        vec .= ∇
        ldiv!(hess, vec)  # `hess \ ∇`, result in `vec`
        @. δx = -vec

        # choosing step in `δx`
        αmax = constrain_step(x, δx)
        α = Downhill.strong_backtracking!(fdf, x, δx;
            α=1.0,  # the function takes care of initial `α` and `αmax`
            αmax=αmax,
            σ=0.9,
        )
        fcalls = 0  # by current implementation they are unknown :c
        totfcalls += fcalls

        if α ≤ 0
            @error "newton: linesearch failed" α
            return NewtonResult(false, x, i, totfcalls)
        end

        # update argument `x`, function value and gradient in `x`
        xpre .= x
        ypre = y

        @. x += α * δx
        y = f(x)
        ∇ = ∇f!(∇, x)

        @debug "newton" i repr(δx) norm_step=norm(δx, 2) α αmax repr(x) y fcalls norm_gprev=norm(∇f!(similar(x), x - α*δx)) norm_gnew=norm(∇) prod(diag(hess.U))

        # check for convergence
        convcond(x, xpre, y, ypre, ∇) && return NewtonResult(true, x, i, totfcalls)
    end
    return NewtonResult(false, x, maxiter, totfcalls)
end

function __convergence_closure(
    state1::AbstractVTFlashState,
    mix::BrusilovskyEoSMixture,
    nmolb::AbstractVector,
    volumeb::Real,
    RT::Real;
    chemtol::Real,
    presstol::Real,
    buf::BrusilovskyThermoBuffer=thermo_buffer(mix),
)
    state1x = value(state1)

    g1 = Vector{Float64}(undef, ncomponents(mix) + 1)
    g2 = similar(g1)
    diff = similar(g1)
    nmol2 = similar(g1, Float64, ncomponents(mix))

    function convcond(x::V, xpre::V, y::T, ypre::T, g::V) where {T<:Real,V<:AbstractVector}
        state1x .= x
        nmol1, vol1 = nmolvol(state1, nmolb, volumeb)
        g1 = CubicEoS.nvtgradient!(g1, mix, nmol1, vol1, RT; buf=buf)

        @. nmol2 = nmolb - nmol1
        vol2 = volumeb - vol1
        g2 = CubicEoS.nvtgradient!(g2, mix, nmol2, vol2, RT; buf=buf)

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
Construct functions of proper signature for [`newton`](@ref) algorithm.
The function reuses `helmholtz_diff!` and `helmholtz_diff_grad!` which are
same as for BFGS-version of vtflash.

(!) BFGS' closures must be constructed from the same thermodynamical properties.
"""
function __vt_flash_newton_closures(
    state1::AbstractVTFlashState,
    helmdiff_qnewton!::Function,
    helmgrad_qnewton!::Function,
    mix::BrusilovskyEoSMixture,
    nmolb::AbstractVector,
    volumeb::Real,
    RT::Real,
)
    auxv = Vector{Float64}(undef, length(nmolb) + 1)
    state1x = value(state1)
    buf = thermo_buffer(mix)

    function clsr_func(x::AbstractVector)
        state1x .= x
        Δa, _ = helmdiff_qnewton!(x, auxv)
        return Δa
    end

    function clsr_hessian!(hess::AbstractMatrix, x::AbstractVector)
        state1x .= x
        hess = hessian!(hess, state1, mix, nmolb, volumeb, RT; buf=buf)
        return hess
    end

    return clsr_func, helmgrad_qnewton!, clsr_hessian!
end

function vt_flash_newton!(
    unstable_state::AbstractVTFlashState,
    mix::BrusilovskyEoSMixture,
    nmol::AbstractVector,
    volume::Real,
    RT::Real;
    chemtol::Real,
    presstol::Real,
    maxiter::Int,
)
    state = unstable_state

    # bfgs closures, just to reuse later
    constrain_step, helmgrad!, helmdiff! = __vt_flash_optim_closures(
        state, mix, nmol, volume, RT
    )

    # closures for newton algorithm
    newton_helmdiff, newton_helmgrad!, newton_helmhess! =
        __vt_flash_newton_closures(
            state,
            helmdiff!,
            helmgrad!,
            mix,
            nmol,
            volume,
            RT,
        )

    convcond = __convergence_closure(state, mix, nmol, volume, RT;
        chemtol=chemtol,
        presstol=presstol,
    )

    statex = value(state)

    optimresult = newton(
        newton_helmdiff,
        newton_helmgrad!,
        newton_helmhess!,
        statex;
        constrain_step=constrain_step,
        gtol=NaN,
        maxiter=maxiter,
        convcond=convcond,
    )
    statex .= optimresult.argument

    return __vt_flash_two_phase_result(state, mix, nmol, volume, RT, optimresult)
end

"""
    vt_flash_newton(mix, nmol, volume, RT, StateVariables[; tol, chemtol, presstol, maxiter])

Find VT-equilibrium of `mix`, at given `nmol`, `volume` and thermal energy `RT`
using Newton's minimization.

For the arguments see [`vt_flash`](@ref).

Return [`VTFlashResult`](@ref).
"""
function vt_flash_newton(
    mix::BrusilovskyEoSMixture{T},
    nmol::AbstractVector,
    volume::Real,
    RT::Real,
    StateVariables::Type{<:AbstractVTFlashState};
    tol::Real=1024*eps(T),
    chemtol::Real=tol,
    presstol::Real=tol,
    maxiter::Int=100,
) where {T}
    singlephase, stability_tries = vt_stability(mix, nmol, volume, RT)

    if singlephase
        concentration = nmol ./ volume
        saturation = 1
        state = StateVariables(concentration, saturation, nmol, volume)
        return __vt_flash_single_phase_result(state, mix, nmol, volume, RT)
    end

    concentration = __vt_flash_init_conc_choose(stability_tries)
    saturation = __find_saturation_negative_helmdiff(mix, nmol, volume, RT, concentration;
        maxsaturation=0.25,
        maxiter=50,
        scale=0.5,
        helmdifftol=-1e-7/RT,
    )

    if isnan(saturation)
        @error "VTFlash: Initial state was not found!" mixture=mix nmol=repr(nmol) volume=volume RT=RT
        error("VTFlash: Initial state was not found!")
    end

    state = StateVariables(concentration, saturation, nmol, volume)

    return vt_flash_newton!(state, mix, nmol, volume, RT;
        chemtol=chemtol,
        presstol=presstol,
        maxiter=maxiter
    )
end
