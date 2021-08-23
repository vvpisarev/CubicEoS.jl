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

    function VTFlashResult{T}(converged, singlephase, RT, nmol_1, V_1, nmol_2, V_2) where {T}
        return new{T}(converged, singlephase, RT, copy(nmol_1), V_1, copy(nmol_2), V_2)
    end
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
            Dmin = state.energy_density
        end
    end
    return vt_stab_tries[index_min].concentration
end

"""
Calculates pressure gradient for Brusilovsky EoS at point (N₁,..., Nₙ, V).
∇P = [∂P/∂Nᵢ..., ∂P/∂V], where i = 1,...,`ncomponents(mix)`.
"""
function __vt_flash_pressure_gradient!(
    ∇P::AbstractVector{T},
    mix::BrusilovskyEoSMixture{T},
    nmol::AbstractVector,
    volume::Real,
    RT::Real,
) where {T}
    # I did not implement this function in src/basic_thermo.jl
    # because the gradient here does not include ∂P/∂T derivative.
    # Maybe, it should be implemented later in src/basic_thermo.jl.

    A, B, C, D, aij = eos_parameters(mix, nmol, RT)

    # hell arithmetics
    # does compiler smart to detect this as constants
    # if plain operation were put in ∂P/∂Nᵢ for-cycle explicitly?
    V = volume  # alias
    VmB⁻¹ = 1 / (V - B)
    ΣnmolbyVmB² = sum(nmol) * VmB⁻¹^2
    DmC = D - C
    VpC⁻¹ = 1 / (V + C)
    VpC⁻² = VpC⁻¹^2
    VpD⁻¹ = 1 / (V + D)
    VpD⁻² = VpD⁻¹^2
    AbyDmC = A / DmC
    VpC⁻¹mVpD⁻¹byDmC² = (VpC⁻¹ - VpD⁻¹) / DmC^2

    # ∂P/∂Nᵢ part
    for (i, substance) in enumerate(components(mix))
        bᵢ, cᵢ, dᵢ = substance.b, substance.c, substance.d
        ∂ᵢA = 2 * dot(nmol, @view aij[:,i])  # ∂A/∂Nᵢ

        ∇P[i] = RT * (VmB⁻¹ + bᵢ * ΣnmolbyVmB²)
                - (
                    (∂ᵢA * DmC - A * (dᵢ - cᵢ)) * VpC⁻¹mVpD⁻¹byDmC²
                    + AbyDmC * (-cᵢ * VpC⁻² + dᵢ * VpD⁻²)
                )
    end
    ∇P[end] = - RT * ΣnmolbyVmB² + AbyDmC * (VpC⁻² - VpD⁻²)
    return nothing
end

"""
Calculates hessian for VTFlash from `state` and base `nmol`, `volume`.
The `state` must be [N₁'/N₁, ..., Nₙ'/Nₙ, V'/V] vector.
"""
function __vt_flash_hessian!(
    hess::AbstractMatrix{T},
    state::AbstractVector{T},
    mix::BrusilovskyEoSMixture{T},
    nmol::AbstractVector,  # I-state
    volume::Real,          # I-state
    RT::Real,
) where {T}
    # TODO: make hessian symmetric for optimization

    #     |    |   |   M   size
    #     | 𝔹  | ℂ |   --------
    # ℍ = |    |   |   𝔹   n×n
    #     |----|---|   ℂ   n×1
    #     | ℂᵀ | 𝔻 |   𝔻   1×1

    #                 [ ∂lnΦᵢ           ∂lnΦᵢ           ]
    # 𝔹ᵢⱼ = -RT Nᵢ Nⱼ [ -----(N', V') + -----(N'', V'') ]
    #                 [  ∂Nⱼ             ∂Nⱼ            ]
    N₁ = nmol .* state[1:end-1]
    V₁ = volume * state[end]

    𝔹 = @view hess[1:end-1, 1:end-1]
    ∇P = similar(state)  # (n + 1) size
    ∇P⁻ = @view ∇P[1:end-1]  # n size
    # ∇P⁻ used as buffer
    log_c_activity_wj!(∇P⁻, 𝔹, mix, N₁, V₁, RT)  # 𝔹 = jacobian'

    N₂ = nmol .- N₁
    V₂ = volume - V₁
    jacobian₂ = Matrix{T}(undef, size(𝔹))
    # ∇P⁻ used as buffer
    log_c_activity_wj!(∇P⁻, jacobian₂, mix, N₂, V₂, RT)

    𝔹 .+= jacobian₂  # 𝔹 = jacobian' + jacobian''
    𝔹 .*= RT * (nmol * nmol')  # final 𝔹, the minus missed cuz of ln Φᵢ = -ln Cₐᵢ

    #            [ ∂P             ∂P             ]
    # ℂᵢ = -V Nᵢ [ --- (N', V') + --- (N'', V'') ]
    #            [ ∂Nᵢ            ∂Nᵢ            ]
    #
    #         [ ∂P            ∂P            ]
    # 𝔻 = -V² [ -- (N', V') + -- (N'', V'') ]
    #         [ ∂V            ∂V            ]
    __vt_flash_pressure_gradient!(∇P, mix, N₁, V₁, RT)
    ℂ = @view hess[1:end-1, end]
    ℂ .= @view ∇P[1:end-1]  # ℂ = (∂P/∂Nᵢ)'
    𝔻 = ∇P[end]  # 𝔻 = (∂P/∂V)'

    __vt_flash_pressure_gradient!(∇P, mix, N₂, V₂, RT)
    ℂ .+= @view ∇P[1:end-1]  # ℂ = ∇P' + ∇P''
    ℂ .*= -volume .* nmol  # final ℂ
    # seems like can be replaced hess[end, :] .= ℂ
    # but version below seems tracking math for me
    hess[[end], 1:end-1] .= ℂ'  # ℂᵀ part of hessian

    𝔻 += ∇P[end]  # 𝔻 = (∂P/∂V)' + (∂P/∂V)''
    𝔻 *= -volume^2  # final 𝔻
    hess[end, end] = 𝔻
    return nothing
end

function vt_flash_closures(
    mix::BrusilovskyEoSMixture{T},
    nmol::AbstractVector,
    volume::Real,
    RT::Real,
) where {T}
    N₁ = Vector{T}(undef, ncomponents(mix))
    N₂ = Vector{T}(undef, ncomponents(mix))
    log_cₐ₁ = Vector{T}(undef, ncomponents(mix))
    log_cₐ₂ = Vector{T}(undef, ncomponents(mix))

    # calculates once
    Pbase = pressure(mix, nmol, volume, RT)
    log_cₐ_base = Vector{T}(undef, ncomponents(mix))
    log_c_activity!(log_cₐ_base, mix, nmol, volume, RT)

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
        if dot(dir, covolumes_b̃) > 0
            if 0 < αm_covolume < αm
                αm = αm_covolume
            end
        else
            if αm < αm_covolume
                @warn "Covolume constrain not meet others"
            end
        end
        if αm == Inf
            error("VTFlash: constrain_step. Step was not found.")
        end
        return 0.9 * αm
    end

    function helmholtz_diff_grad!(state::AbstractVector{T}, grad::AbstractVector{T})
        _, V₁, V₂ = transform(state)
        log_c_activity!(log_cₐ₁, mix, N₁, V₁, RT)
        log_c_activity!(log_cₐ₂, mix, N₂, V₂, RT)

        @inbounds for i in 1:length(state)-1
            Δμ = -RT * (log((N₂[i]/V₂) / (N₁[i]/V₁)) - (log_cₐ₁[i] - log_cₐ₂[i]))
            grad[i] = nmol[i] * Δμ
        end
        P₁ = pressure(mix, N₁, V₁, RT)
        P₂ = pressure(mix, N₂, V₂, RT)
        grad[end] = volume * (-P₁ + P₂)
        return grad
    end
    function helmholtz_diff!(state::AbstractVector{T}, grad::AbstractVector{T})
        _, V₁, V₂ = transform(state)

        log_c_activity!(log_cₐ₂, mix, N₂, V₂, RT)

        "Σ Nᵢ (μᵢ - μ₂ᵢ)"
        Ndotμ₂ = zero(T)
        @inbounds for i in 1:length(state)-1
            # μ base - μ₂
            Δμ = -RT * (log((N₂[i]/V₂)/(nmol[i]/volume)) - (log_cₐ_base[i] - log_cₐ₂[i]))
            Ndotμ₂ += nmol[i] * Δμ
        end

        P₂ = pressure(mix, N₂, V₂, RT)
        helmholtz_diff_grad!(state, grad)  # overwrites gradient `grad`
        ΔA = dot(grad, state) + (Pbase - P₂) * volume - Ndotμ₂
        @debug "helmholtz_diff!" state=repr(state) ΔA grad=repr(grad) norm(grad, 2)
        return ΔA, grad
    end
    return constrain_step, helmholtz_diff_grad!, helmholtz_diff!
end

function __vt_flash_initial_state!(
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
    @debug "Initial state search" start_scale=scale sat₁max
    for i in 1:steps
        # upd `state`
        sat = sat₁max * scale
        state[1:end-1] .= conc₁ * (sat * volume) ./ nmol
        state[end] = sat
        # TODO: check if state feasible

        # calc helmholtz energy
        @debug "Initial state search" i state=repr(state) scale
        try
            ΔA, _ = helmholtz_diff!(state, vec)
            # check convergence
            @debug "Initial state search: ΔA calculated" i ΔA helmholtz_thresh
            if ΔA < helmholtz_thresh
                return true
            end
        catch e
            @warn "VTFlash: initial state search" sat e
        end

        # update `scale`
        scale *= step_scale
    end
    return false
end

"""
Extracts vt-state from `optresult` (DescentMethods obj).
Sorts variables into gas and liquid.
Returns corresponding `VTFlashResult`.
"""
function __vt_flash_two_phase_result(
    mix::BrusilovskyEoSMixture{T},
    nmol::AbstractVector{T},
    volume::Real,
    RT::Real,
    optresult,
) where {T}
    # □₁ for gas, □₂ for liquid
    state = optresult.argument
    nmol₁ = nmol .* @view state[1:end-1]
    V₁ = volume * state[end]
    nmol₂ = nmol .- nmol₁
    V₂ = volume - V₁

    P₁ = pressure(mix, nmol₁, V₁, RT)  # they should be equal
    P₂ = pressure(mix, nmol₂, V₂, RT)

    Z₁ = P₁ * V₁ / (sum(nmol₁) * RT)  # seems can be reduced to Vᵢ / sum(nmolᵢ)
    Z₂ = P₂ * V₂ / (sum(nmol₂) * RT)

    if Z₂ > Z₁  # □₂ is gas state, need exchange
        P₁, P₂ = P₂, P₁
        Z₁, Z₂ = Z₂, Z₁
        V₁, V₂ = V₂, V₁

        for i in eachindex(nmol₁, nmol₂)
            nmol₁[i], nmol₂[i] = nmol₂[i], nmol₁[i]
        end
    end

    return VTFlashResult{T}(
            converged=optresult.converged,
            singlephase=false,
            RT=RT,
            nmol_1=nmol₁,
            V_1=V₁,
            nmol_2=nmol₂,
            V_2=V₂
    )
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

    # initial hessian
    hessian = Matrix{T}(undef, (length(state), length(state)))
    __vt_flash_hessian!(hessian, state, mix, nmol, volume, RT)
    @debug "VTFlash: initial hessian found" isposdef(hessian)

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

    return __vt_flash_two_phase_result(mix, nmol, volume, RT, result)
end
