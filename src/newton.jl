using PositiveFactorizations
using LinearAlgebra


struct NewtonResult{T}
    converged::Bool
    argument::Vector{T}
    NewtonResult(conv, x) = new{Float64}(conv, copy(x))
end

"""
`∇f(x) -> gradient` return gradient of `f` at `x`.
`H -> hessian` return hessian of `f` at `x`.
`x` -> initial `x`.
"""
function newton(∇f, H, x;
    ∇atol=1e-6,
    maxiter=100,
    constrain_step=(x, δ)->1,
)
    for i in 1:maxiter
        H̃ = cholesky(Positive, H(x))
        ∇ = ∇f(x)
        δx = - (H̃ \ ∇)  # brackets because there is no unary - for ::Cholesky{...}
        α = constrain_step(x, δx)

        @debug "newton" i repr(δx) repr(α) repr(x) norm(∇f(x)) isposdef(H̃)

        x += min(1, α) * δx
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

    # gradient function
    gradient = Vector{T}(undef, length(state))
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
    result = newton(helmholtz_diff_gradient, helmholtz_diff_hessian, state;
        constrain_step=constrain_step,
        ∇atol=1e-3,
    )

    return __vt_flash_two_phase_result(mix, nmol, volume, RT, result)
end
