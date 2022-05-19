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

Similar to [`vt_flash`](@ref), find thermodynamical equilibrium of `mix`
at given `nmol`, `volume` and thermal energy `RT`, but uses *Newton's optimization* in
phase-split stage.

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
    stabconverged, singlephase, stability_tries = vt_stability(mix, nmol, volume, RT)

    !stabconverged && error("VT-Stability does not converged")

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
