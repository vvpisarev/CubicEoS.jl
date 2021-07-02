#=
VT-flash algorithm

state assumed to be [molar fractions of phase 1, saturation of phase 1]
=#

struct VTFlashResult{T}
    converged::Bool
    singlephase::Bool
    RT::T
    nmol_1::Vector{T}
    V_1::T
    nmol_2::Vector{T}
    V_2::T
end

VTFlashResult{T}(; converged, singlephase, RT, nmol_1, V_1, nmol_2, V_2) where {T} =
VTFlashResult{T}(converged, singlephase, RT, nmol_1, V_1, nmol_2, V_2)

function vt_flash_closures(
    mix::BrusilovskyEoSMixture{T},
    nmol::AbstractVector,
    volume::Real,
    RT::Real,
) where {T}
    N₁ = Vector{T}(undef, ncomponents(mix))
    N₂ = Vector{T}(undef, ncomponents(mix))
    log_Φ₁ = Vector{T}(undef, ncomponents(mix))
    log_Φ₂ = Vector{T}(undef, ncomponents(mix))

    covolumes_b̃ = [(c.b for c in components(mix))..., -volume/sum(nmol)]

    # calculates once
    Σnmol = sum(nmol)
    Pbase = pressure(mix, nmol, volume, RT)
    log_Φbase = Vector{T}(undef, ncomponents(mix))
    log_c_activity!(log_Φbase, mix, nmol, volume, RT)

    function constrain_step(state::AbstractVector{T}, dir::AbstractVector{T})
        αm = convert(T, Inf)
        # positiveness constrain
        @inbounds for i in eachindex(state)
            α = (1 - state[i]) / dir[i]
            if 0 < α < αm
                αm = α
            elseif α < 0
                @warn "constrain not meet for i = $i" state[i] dir[i] α
            end
        end
        # covolume constrain
        αm_covolume = - dot(state, covolumes_b̃) / dot(dir, covolumes_b̃)
        if 0 < αm_covolume < αm
            αm = αm_covolume
        end
        return αm
    end

    function helmholtz_diff_grad!(state::AbstractVector{T}, grad_::AbstractVector{T})
        molfrac₁ = @view state[1:end-1]
        sat₁ = state[end]

        # log_Φ₁
        N₁ .= Σnmol .* molfrac₁
        V₁ = volume * sat₁
        log_c_activity!(log_Φ₁, mix, N₁, V₁, RT)
        # log_Φ₂
        N₂ .= nmol .- N₁
        V₂ = volume - V₁
        log_c_activity!(log_Φ₂, mix, N₂, V₂, RT)

        @inbounds for i in eachindex(molfrac₁)
            χ₁, S₁ = molfrac₁[i], sat₁
            χ₂, S₂ = 1 - χ₁, 1 - sat₁
            # (χ₂/S₂)/(χ₁/S₁) ≡ (N₂/V₂) / (N₁/V₁)
            Δμ = -RT * (log((χ₂/S₂)/(χ₁/S₁)) + (log_Φ₁[i] - log_Φ₂[i]))
            grad_[i] = Δμ
        end
        P₁ = pressure(mix, N₁, V₁, RT)
        P₂ = pressure(mix, N₂, V₂, RT)
        grad_[end] = -P₁ + P₂
        return grad_
    end
    function helmholtz_diff!(state::AbstractVector{T}, grad_::AbstractVector{T})
        helmholtz_diff_grad!(state, grad_)  # overwrites gradient `grad_`

        N₂ .= nmol .- (sum(nmol) .* @view state[1:end-1])
        V₂ = volume * (1 - state[end])
        log_c_activity!(log_Φ₂, mix, N₂, V₂, RT)

        "Σ Nᵢ (μᵢ - μ₂ᵢ)"
        Ndotμ₂ = zero(T)
        @inbounds for i in 1:length(state)-1
            # μ base - μ₂
            χ, S = nmol[i] / Σnmol, 1
            χ₂, S₂ = 1 .- state[i], 1 - state[end]
            Δμ = -RT * (log((χ₂/S₂)/(χ/S)) + (log_Φbase[i] - log_Φ₂[i]))
            Ndotμ₂ += nmol[i] * Δμ
        end

        P₂ = pressure(mix, N₂, V₂, RT)
        ΔA = dot(grad_, state) + (Pbase - P₂) * volume - Ndotμ₂
        return ΔA, grad_
    end
    return constrain_step, helmholtz_diff_grad!, helmholtz_diff!
end

function vt_flash_initial_state!(
    state::AbstractVector{T},
    conc₁::AbstractVector{T},
    helmholtz_diff!::Function;
    sat₁max::Real=T(0.5),
    steps::Int=20,
    step_scale::Real=T(0.5),
    helmholtz_thresh::Real=T(-1e-5),  # must be negative value
) where {T}
    for i in 1:length(state)-1
        state[i] = conc₁[i] / sum(conc₁)
    end
    state[end] = sat₁max

    vec = similar(state)  # buffer vector for gradient
    scale = 1
    for i in 1:steps
        # upd `state`
        state[end] *= scale
        # TODO: check if state feasible

        # calc helmholtz energy
        ΔA, _ = helmholtz_diff!(state, vec)
        # check convergence
        @debug i state=repr(state) scale ΔA helmholtz_thresh
        if ΔA < helmholtz_thresh
            return state
        end

        # update `scale`
        scale *= step_scale
    end
    error("vt_flash_initial_state: state not found")
    return nothing
end

function vt_flash(
    mix::BrusilovskyEoSMixture{T},
    nmol::AbstractVector,
    volume::Real,
    RT::Real,
) where {T}
    # run vt-stability to find out whether a state single phase or not
    singlephase, η_base = vt_stability(mix, nmol, volume, RT)
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
    vt_flash_initial_state!(
        state,
        η_base,
        helmholtz_diff!;
        sat₁max=0.9,
        steps=20,
        step_scale=0.5,
        helmholtz_thresh=-1e-5,
    )

    # initial hessian
    hessian = 1e-7 * ones((length(state), length(state)))

    # run optimizer
    optmethod = DescentMethods.CholBFGS(state)
    DescentMethods.reset!(optmethod, state, hessian)
    result = DescentMethods.optimize!(
        optmethod,
        helmholtz_diff!,
        state,
        gtol=1e-3,
        maxiter=1000,
        constrain_step=constrain_step,
        reset=false,
    )
    # TODO: check convergence of BFGS

    # TODO: sort phases into foo_1 for gas, foo_2 for liquid
    state .= result.x
    nmol₁ = nmol .* @view state[1:end-1]
    V₁ = volume * state[end]
    nmol₂ = nmol .- nmol₁
    V₂ = volume - V₁

    return VTFlashResult{T}(
            converged=result.converged,
            singlephase=false,
            RT=RT,
            nmol_1=nmol₁,
            V_1=V₁,
            nmol_2=nmol₂,
            V_2=V₂
    )
end
