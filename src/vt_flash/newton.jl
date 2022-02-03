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
    x::AbstractVector;
    gtol::Real=1e-6,
    maxiter::Integer=200,
    constrain_step::Function=(x, δ)->Inf,
)
    x = convert(Vector{Float64}, x)
    ∇ = similar(x)
    δx = similar(x)
    hess_full = Matrix{Float64}(undef, size(∇, 1), size(∇, 1))
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
    if norm(∇, 2) ≤ gtol
        return NewtonResult(true, x, 0, totfcalls)
    end

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
        if norm(∇, 2) ≤ gtol
            return NewtonResult(true, x, i, totfcalls)
        end
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
    helmholtz_diff!::Function,
    helmholtz_diff_grad!::Function,
    mix::BrusilovskyEoSMixture{T},
    nmol::AbstractVector{T},
    volume::T,
    RT::T,
) where {T<:Real}
    vecnc₊ = Vector{T}(undef, size(nmol, 1) + 1)
    thermo_buf = thermo_buffer(mix)
    hessian_buf = HessianBuffer(mix)

    # below `x` is mixtures' state
    # helmholtz
    function func(x::AbstractVector{T})
        ΔA, _ = helmholtz_diff!(x, vecnc₊)
        return ΔA
    end

    # helmholtz gradient
    function gradient!(gradient_::AbstractVector{T}, x::AbstractVector{T})
        helmholtz_diff_grad!(x, gradient_)  # yes, base func uses this signature
        return gradient_
    end

    # helmholtz hessian
    function hessian!(hessian::AbstractMatrix{T}, x::AbstractVector{T})
        __vt_flash_hessian!(hessian, x, mix, nmol, volume, RT; buf=hessian_buf)
        return hessian
    end

    return func, gradient!, hessian!
end

function vt_flash_newton(
    mix::BrusilovskyEoSMixture{T},
    nmol::AbstractVector,
    volume::Real,
    RT::Real,
    unstable_state::AbstractVector,
) where {T}
    state = copy(unstable_state)

    # bfgs closures, just to reuse later
    constrain_step, helmholtz_diff_grad!, helmholtz_diff! = vt_flash_closures(
        mix, nmol, volume, RT
    )
    # closures for newton algorithm
    newton_helmholtz_diff, newton_helmholtz_diff_grad!, newton_helmholtz_diff_hessian! =
        __vt_flash_newton_closures(
            helmholtz_diff!,
            helmholtz_diff_grad!,
            mix,
            nmol,
            volume,
            RT,
        )

    # run optimizer
    result = newton(
        newton_helmholtz_diff,
        newton_helmholtz_diff_grad!,
        newton_helmholtz_diff_hessian!,
        state;
        constrain_step=constrain_step,
        gtol=1e-3,
        maxiter=200,
    )

    return __vt_flash_two_phase_result(mix, nmol, volume, RT, result)
end

"""
    vt_flash_newton(mix, nmol, volume, RT)

Find VT-equilibrium of `mix`, at given `nmol`, `volume` and thermal energy `RT`
using Newton's minimization. Return [`VTFlashResult`](@ref).
"""
function vt_flash_newton(
    mix::BrusilovskyEoSMixture{T},
    nmol::AbstractVector,
    volume::Real,
    RT::Real,
) where {T}
    # run vt-stability to find out whether a state single phase or not
    singlephase, vt_stab_tries = vt_stability(mix, nmol, volume, RT)
    @debug "VTFlash: VTStability result" singlephase

    if singlephase
        return VTFlashResult{T}(
            converged=true,
            singlephase=true,
            RT=RT,
            nmol_1=nmol,
            V_1=volume,
            nmol_2=similar(nmol),
            V_2=0
        )
    end

    # two-phase state case
    init_found, state = __vt_flash_initial_state(
        mix, nmol, volume, RT, vt_stab_tries
    )

    @debug "VTFlash: initial state search result" found=init_found state=repr(state) ΔA=helmholtz_diff!(state, similar(state))

    if !init_found
        @error "VTFlash: Initial state was not found!" mixture=mix nmol=repr(nmol) volume=volume RT=RT
        error("VTFlash: Initial state was not found!")
    end

    return vt_flash_newton(mix, nmol, volume, RT, state)
end
