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
    newton(f, grad!, hess!, x₀[; gtol, maxiter, constrain_step])

Find constrained minimum of `f(x)::Real` using Newton's solver with line search
and hessian modified by Cholesky factorization.
`grad!(g, x)` must update and return `g`radient of `f` at `x`.
`hess!(H, x)` must update and return `H`essian matrix of `f` at `x`.

The optimization starts at `x₀` performing no more than `maxiter` Newton's steps.
The steps finish when norm of gradient is ≤ `gtol`.
New point `xnew <- x + αmax * d` is constrained in Newton's `d`irection by
`constrain_step(x, d) -> α::Real`, `αmax = min(1.0, α)`.

Return [`NewtonResult`](@ref) object.
"""
function newton(
    f::Function,
    ∇f!::Function,
    H!::Function,
    x_::AbstractVector;
    gtol::Real=1e-6,
    maxiter::Integer=200,
    constrain_step::Function=(x, δ)->Inf,
)
    x = copy(x_)  # copy x for internal use

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

    # check gradient norm before iterations
    ∇ = ∇f!(∇, x)
    norm(∇, 2) ≤ gtol && return NewtonResult(true, x, 0, totfcalls)

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
        norm(∇, 2) ≤ gtol && return NewtonResult(true, x, i, totfcalls)
    end
    return NewtonResult(false, x, maxiter, totfcalls)
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
    gtol::Real,
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

    statex = value(state)

    optimresult = newton(
        newton_helmdiff,
        newton_helmgrad!,
        newton_helmhess!,
        statex;
        constrain_step=constrain_step,
        gtol=gtol,
        maxiter=maxiter,
    )
    statex .= optimresult.argument

    return __vt_flash_two_phase_result(state, mix, nmol, volume, RT, optimresult)
end

"""
    vt_flash_newton(mix, nmol, volume, RT)

Find VT-equilibrium of `mix`, at given `nmol`, `volume` and thermal energy `RT`
using Newton's minimization. Return [`VTFlashResult`](@ref).
"""
function vt_flash_newton(
    mix::BrusilovskyEoSMixture,
    nmol::AbstractVector,
    volume::Real,
    RT::Real,
    StateVariables::Type{<:AbstractVTFlashState};
    gtol::Real=1e-3/RT,
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

    return vt_flash_newton!(state, mix, nmol, volume, RT; gtol=gtol, maxiter=maxiter)
end