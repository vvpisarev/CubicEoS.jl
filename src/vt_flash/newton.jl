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
    newton(f, grad!, hess!, x₀[; gtol=1e-6, maxiter=200, constrain_step, convcond])

Find (un)constrained minimum of `f(x)::Real` using Newton's solver with line search
and hessian modified by Cholesky factorization.

The optimization starts at `x₀` performing no more than `maxiter` Newton's steps.

New point `xnew ← x + αmax * d` is constrained in Newton's `d`irection by
`constrain_step(x, d) -> α::Real`.
By default, the constrains are disabled.

The *end of optimization* is controlled by `gtol` and `convcond`.
When `convcond` is omitted, the optimization ends when 2-norm of gradient is ≤ `gtol`.
When `convcond(x, g) -> Bool` is provided, it's called before checking of gradient's norm.

Return [`NewtonResult`](@ref) object.

# Arguments

- `grad!(g, x) -> g`: must update and return `g`radient of `f` at `x`;
- `hess!(H, x) -> H`: must update and return `H`essian matrix of `f` at `x`;

# Optional arguments

- `gtol::Real=1e-6`: 2-norm of gradient when iterations are considered converged,
                     if `NaN`, the norm isn't checked;
- `maxiter::Integer=200`: maximum allowed Newton's steps;
- `constrain_step::Function=(x, d)->Inf`:
        given an argument `x` and `d`irection, returns maximum allowed scalar `α` in form of
        `xnew = x + α*d`. Default is unconstrained optimization.
- `convcond::Function=(xnew, gnew)->false`:
        if returns `true` for given new point and gradient, the optimization is stopped.
        Disabled by default.
"""
function newton(
    f::Function,
    ∇f!::Function,
    H!::Function,
    x_::AbstractVector;
    gtol::Real=1e-6,
    maxiter::Integer=200,
    constrain_step::Function=(x, d)->Inf,
    convcond::Function=(xnew, gnew)->false,
)
    x = float.(x_)  # copy x for internal use

    ∇ = similar(x)
    δx = similar(x)
    hess_full = Matrix{eltype(x)}(undef, size(∇, 1), size(∇, 1))
    vec = similar(x)

    fval = f(x)
    totfcalls = 1

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

    # check for convergence in initial point
    ∇ = ∇f!(∇, x)
    convcond(x, ∇) && return NewtonResult(true, x, 0, totfcalls)
    !isnan(gtol) && norm(∇, 2) ≤ gtol && return NewtonResult(true, x, 0, totfcalls)

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

        # previous code, which uses "simple" backtracking
        # αmax = min(1.0, constrain_step(x, δx))
        # α, fval, fcalls = backtracking_line_search(f, x, δx, fval; α=αmax, p=0.5, buf=vec)

        totfcalls += fcalls

        if α ≤ 0
            @error "newton: linesearch failed" α
            return NewtonResult(false, x, i, totfcalls)
        end

        # update argument `x`, function value and gradient in `x`
        @. x += α * δx
        fval = f(x)
        ∇ = ∇f!(∇, x)

        @debug "newton" i repr(δx) norm_step=norm(δx, 2) α αmax repr(x) fval fcalls norm_gprev=norm(∇f!(similar(x), x - α*δx)) norm_gnew=norm(∇) prod(diag(hess.U))

        # check for convergence
        convcond(x, ∇) && return NewtonResult(true, x, i, totfcalls)
        !isnan(gtol) && norm(∇, 2) ≤ gtol && return NewtonResult(true, x, i, totfcalls)
    end
    return NewtonResult(false, x, maxiter, totfcalls)
end

function __convergence_closure(
    state1::AbstractVTFlashState,
    mix::BrusilovskyEoSMixture,
    nmolb::AbstractVector,
    volumeb::Real,
    RT::Real;
    chemrtol::Real,
    pressrtol::Real,
    buf::BrusilovskyThermoBuffer=thermo_buffer(mix),
)
    state1x = value(state1)

    g1 = Vector{Float64}(undef, ncomponents(mix) + 1)
    g2 = similar(g1)
    diff = similar(g1)
    nmol2 = similar(g1, Float64, ncomponents(mix))

    function convcond(x::AbstractVector, g::AbstractVector)
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

        nrmchem = let chem1 = (@view g1[1:end-1]), diffchem = (@view diff[1:end-1])
            abs.(diffchem) ./ abs.(chem1) |> maximum
        end
        nrmpress = abs(diff[end] / g1[end])

        return nrmchem ≤ chemrtol && nrmpress ≤ pressrtol
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
    chemrtol::Real,
    pressrtol::Real,
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
        chemrtol=chemrtol,
        pressrtol=pressrtol,
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
    vt_flash_newton(mix, nmol, volume, RT, StateVariables[; chemrtol, pressrtol, maxiter])

Find VT-equilibrium of `mix`, at given `nmol`, `volume` and thermal energy `RT`
using Newton's minimization.

Return [`VTFlashResult`](@ref).
"""
function vt_flash_newton(
    mix::BrusilovskyEoSMixture,
    nmol::AbstractVector,
    volume::Real,
    RT::Real,
    StateVariables::Type{<:AbstractVTFlashState};
    chemrtol::Real=1e3*eps(),
    pressrtol::Real=1e3*eps(),
    maxiter::Int=100,
)
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
        chemrtol=chemrtol,
        pressrtol=pressrtol,
        maxiter=maxiter
    )
end
