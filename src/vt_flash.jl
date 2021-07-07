#=
VT-flash algorithm
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

"Return concentration of state with minimum energy from vt-stability tries."
function __vt_flash_init_conc_choose(
    vt_stab_tries::AbstractVector{VTStabilityResult{T}},
) where {T}
    Dmin = T(Inf)
    index_min = -1
    for (i, state) in enumerate(vt_stab_tries)
        if !state.isstable && state.energy_density < Dmin
            index_min = i
        end
    end
    return vt_stab_tries[index_min].concentration
end

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

    # calculates once
    Σnmol = sum(nmol)
    Pbase = pressure(mix, nmol, volume, RT)
    log_Φbase = Vector{T}(undef, ncomponents(mix))
    log_c_activity!(log_Φbase, mix, nmol, volume, RT)

    "Constant vector for covolume constrain. [Nᵢbᵢ..., -V]"
    covolumes_b̃ = [(c.b for c in components(mix))..., 1]
    covolumes_b̃[1:end-1] .*= nmol
    covolumes_b̃[end] *= -volume

    "Updates `N₁`, `N₂`. Returns `state_tr`, `V₁`, `V₂` from `state`."
    function transform(state::AbstractVector{T})
        Tr = Diagonal([nmol..., volume])
        state_tr = Tr * state
        N₁ .= state_tr[1:end-1]
        N₂ .= nmol .- N₁
        V₁ = state_tr[end]
        V₂ = volume - V₁
        return state_tr, V₁, V₂
    end

    function constrain_step(state::AbstractVector{T}, dir::AbstractVector{T})
        αm = convert(T, Inf)
        # # positiveness constrain (`0 < state[i] + α * dir[i] < 1`)
        @inbounds for i in eachindex(state)
            if dir[i] > 0
                α = (1 - state[i]) / dir[i]
            elseif dir[i] < 0
                α = - state[i] / dir[i]
            end
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
        if αm == Inf
            error("VTFlash: constrain_step. Step was not found.")
        end
        return 0.9 * αm
    end

    function helmholtz_diff_grad!(state::AbstractVector{T}, grad_::AbstractVector{T})
        _, V₁, V₂ = transform(state)
        log_c_activity!(log_Φ₁, mix, N₁, V₁, RT)
        log_c_activity!(log_Φ₂, mix, N₂, V₂, RT)

        @inbounds for i in 1:length(state)-1
            Δμ = -RT * (log((N₂[i]/V₂) / (N₁[i]/V₁)) + (log_Φ₁[i] - log_Φ₂[i]))
            grad_[i] =  nmol[i] * Δμ
        end
        P₁ = pressure(mix, N₁, V₁, RT)
        P₂ = pressure(mix, N₂, V₂, RT)
        grad_[end] = volume * (-P₁ + P₂)
        return grad_
    end
    function helmholtz_diff!(state::AbstractVector{T}, grad_::AbstractVector{T})
        _, V₁, V₂ = transform(state)

        # @debug "helmholtz_diff!" state=repr(state) N=repr(nmol) V=volume N₁=repr(N₁) V₁ N₂=repr(N₂) V₂
        # @debug "helmholtz_diff!: constrains" dot(N₁, covolumes_b̃[1:end-1]./nmol) dot(N₂, covolumes_b̃[1:end-1]./nmol)

        log_c_activity!(log_Φ₂, mix, N₂, V₂, RT)

        "Σ Nᵢ (μᵢ - μ₂ᵢ)"
        Ndotμ₂ = zero(T)
        @inbounds for i in 1:length(state)-1
            # μ base - μ₂
            Δμ = -RT * (log((N₂[i]/V₂)/(nmol[i]/volume)) + (log_Φbase[i] - log_Φ₂[i]))
            Ndotμ₂ += nmol[i] * Δμ
        end

        P₂ = pressure(mix, N₂, V₂, RT)
        helmholtz_diff_grad!(state, grad_)  # overwrites gradient `grad_`
        ΔA = dot(grad_, state) + (Pbase - P₂) * volume - Ndotμ₂
        @debug "helmholtz_diff!" state=repr(state) ΔA grad_ΔA=repr(grad_)
        return ΔA, grad_
    end
    return constrain_step, helmholtz_diff_grad!, helmholtz_diff!
end

function vt_flash_initial_state!(
    state::AbstractVector{T},
    nmol::AbstractVector{T},
    volume::Real,
    conc₁::AbstractVector{T},
    helmholtz_diff!::Function,
    constrain_step::Function;
    sat₁max::Real=T(0.5),
    steps::Int=20,
    step_scale::Real=T(0.5),
    helmholtz_thresh::Real=T(-1e-5),  # must be negative value
) where {T}
    state[1:end-1] .= conc₁ * (sat₁max * volume) ./ nmol
    state[end] = sat₁max

    vec = similar(state)  # buffer vector for gradient
    scale = one(T)
    # scale = 0.9 * constrain_step(state, -1 * ones(length(state)))
    @debug "Initial state search" start_scale=scale sat₁max
    for i in 1:steps
        # upd `state`
        sat = sat₁max * scale
        state[1:end-1] .= conc₁ * (sat * volume) ./ nmol
        state[end] = sat
        # TODO: check if state feasible

        # calc helmholtz energy
        @debug "Initial state search" i state=repr(state) scale
        ΔA, _ = helmholtz_diff!(state, vec)
        # check convergence
        @debug "Initial state search: ΔA calculated" i ΔA helmholtz_thresh
        if ΔA < helmholtz_thresh
            return true
        end

        # update `scale`
        scale *= step_scale
    end
    return false
end

function vt_flash(
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

    init_found = vt_flash_initial_state!(
        state,
        nmol,
        volume,
        η₁test,
        helmholtz_diff!,
        constrain_step;
        sat₁max=0.1,
        steps=20,
        step_scale=0.5,
        helmholtz_thresh=-1e-5,
    )

    @debug "VTFlash: initial state search result" found=init_found state=repr(state) ΔA=helmholtz_diff!(state, similar(state))

    if !init_found
        error("VTFlash: Initial state was not found!")
    end

    # initial hessian
    hessian = 1e-8 * ones((length(state), length(state)))

    # run optimizer
    optmethod = DescentMethods.CholBFGS(state)
    DescentMethods.reset!(optmethod, state, hessian)
    result = DescentMethods.optimize!(
        optmethod,
        helmholtz_diff!,
        state,
        gtol=1e-3,
        maxiter=100,
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
