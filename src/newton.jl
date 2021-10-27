using PositiveFactorizations
using LinearAlgebra

struct Newton{M<:AbstractMatrix}
    hess::M
end

struct NewtonResult{T}
    converged::Bool
    argument::Vector{T}
    iterations::Int
    calls::Int

    NewtonResult(conv, x, iters=-1, calls=-1) = new{Float64}(conv, copy(x), iters, calls)
end

function backtracking_line_search(
    f::Function,
    x::AbstractVector{T},
    d::AbstractVector{T},
    α::T=one(T),
    p::Real=one(T)/2,
    buffer::AbstractVector{T}=similar(x),
) where {T<:Real}
    calls = 1
    y = f(x)
    α = T(α)
    xtry = buffer
    while true
        @. xtry = x + α*d
        calls += 1
        f(xtry) < y && break
        α *= p
        @debug "backtracking_line_search" α p
    end
    return α, calls
end

"""
`∇f!(∇, x) -> update ∇ with gradient and return it
`H!(hess, x) -> hessian` return hessian of `f` at `x`.
`x` -> initial `x`.
"""
function newton(
    f::Function,
    ∇f!::Function,
    H!::Function,
    x::AbstractVector;
    ∇atol::Real=1e-6,
    maxiter::Integer=200,
    constrain_step::Function=(x, δ)->one(Float64),
)
    x = convert(Vector{Float64}, x)
    ∇ = similar(x)
    δx = similar(x)
    hess_full = Matrix{Float64}(undef, size(∇, 1), size(∇, 1))
    vec = similar(x)
    totfcalls = 0

    for i in 1:maxiter
        hess_full = H!(hess_full, x)  # 6 allocs
        hess = DescentMethods.mcholesky!(hess_full)

        ∇ = ∇f!(∇, x)

        vec .= ∇
        ldiv!(hess, vec)  # `hess \ ∇`, result in `vec`
        @. δx = -vec

        αmax = min(0.5, constrain_step(x, δx))
        α, fcalls = backtracking_line_search(f, x, δx, αmax, 0.5, vec)

        @. x += α * δx

        @debug "newton" i repr(δx) norm(δx, 2) α repr(x) f(x) fcalls norm(∇, 2) det(hess)

        totfcalls += fcalls

        if norm(∇, 2) ≤ ∇atol
            return NewtonResult(true, x, i, totfcalls)
        end
    end
    return NewtonResult(false, x, maxiter, totfcalls)
end

function __vt_flash_newton_closures(
    helmholtz_diff!::Function,      # bfgs vtflash funcs
    helmholtz_diff_grad!::Function, #
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
    # create closures for helmoltz energy, its gradient and constrain step
    constrain_step, helmholtz_diff_grad!, helmholtz_diff! = vt_flash_closures(
        mix, nmol, volume, RT
    )

    # find initial vector for optimizer
    state = Vector{T}(undef, ncomponents(mix) + 1)
    η₁test = __vt_flash_init_conc_choose(vt_stab_tries)

    init_found = __vt_flash_initial_state!(
        state, nmol, volume, η₁test, helmholtz_diff!, constrain_step;
        sat₁max=0.25,
        steps=200,
        step_scale=0.5,
        helmholtz_thresh=-1e-7,
    )

    @debug "VTFlash: initial state search result" found=init_found state=repr(state) ΔA=helmholtz_diff!(state, similar(state))

    if !init_found
        @error "VTFlash: Initial state was not found!" mixture=mix nmol=repr(nmol) volume=volume RT=RT
        error("VTFlash: Initial state was not found!")
    end

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

    @debug "VTFlash: initial" isposdef(helmholtz_diff_hessian(state))

    # run optimizer
    result = newton(
        newton_helmholtz_diff,
        newton_helmholtz_diff_grad!,
        newton_helmholtz_diff_hessian!,
        state;
        constrain_step=constrain_step,
        ∇atol=1e-3,
        maxiter=200,
    )

    return __vt_flash_two_phase_result(mix, nmol, volume, RT, result)
end
