using PositiveFactorizations
using LinearAlgebra


struct NewtonResult{T}
    converged::Bool
    argument::Vector{T}
    NewtonResult(conv, x) = new{Float64}(conv, copy(x))
end

function backtracking_line_search(f, ∇f, x, d, α; p=0.5, β=1e-4)
    y, g = f(x), ∇f(x)
    while f(x + α*d) > y
        α *= p
        @debug "backtracking_line_search" α p
    end
    return α
end

"""
`∇f(x) -> gradient` return gradient of `f` at `x`.
`H -> hessian` return hessian of `f` at `x`.
`x` -> initial `x`.
"""
function newton(f, ∇f, H, x;
    ∇atol=1e-6,
    maxiter=1000,
    constrain_step=(x, δ)->1,
)
    for i in 1:maxiter
        hess = copy(H(x))
        if !isposdef(hess)
            H̃ = DescentMethods.mcholesky!(hess)
        else
            H̃ = hess
        end
        ∇ = ∇f(x)
        δx = - (H̃ \ ∇)

        αmax = min(0.5, constrain_step(x, δx))
        α = backtracking_line_search(f, ∇f, x, δx, αmax; p=0.5, β=1e-4)

        x += α * δx

        @debug "newton" i repr(δx) norm(δx, 2) α repr(x) f(x) norm(∇f(x)) isposdef(H̃) det(H̃)

        if norm(∇, 2) ≤ ∇atol
            return NewtonResult(true, x)
        end
    end
    return NewtonResult(false, x)
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

    gradient = Vector{T}(undef, length(state))
    helmholtz_diff = function(x)
        ΔA, _ = helmholtz_diff!(x, gradient)
        return ΔA
    end
    # gradient function
    helmholtz_diff_gradient = function (x)
        helmholtz_diff_grad!(x, gradient)
        return gradient
    end
    # hessian function
    hessian = Matrix{T}(undef, (length(state), length(state)))
    helmholtz_diff_hessian = function (x)
        __vt_flash_hessian!(hessian, x, mix, nmol, volume, RT)
        return hessian
    end

    @debug "VTFlash: initial" isposdef(helmholtz_diff_hessian(state))

    # run optimizer
    result = newton(
        helmholtz_diff,
        helmholtz_diff_gradient,
        helmholtz_diff_hessian,
        state;
        constrain_step=constrain_step,
        ∇atol=1e-3,
    )

    return __vt_flash_two_phase_result(mix, nmol, volume, RT, result)
end
