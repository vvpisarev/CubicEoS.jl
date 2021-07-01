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
    nc = ncomponents(mix)
    N₁ = Vector{T}(undef, ncomponents(mix))
    N₂ = Vector{T}(undef, ncomponents(mix))
    log_Φ₁ = Vector{T}(undef, ncomponents(mix))
    log_Φ₂ = Vector{T}(undef, ncomponents(mix))
    Σnmol = sum(nmol)

    function constrain_step(state, dir) end
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
        grad_[nc+1] = -P₁ + P₂
        return grad_
    end
    function helmholtz_diff!(state::AbstractVector{T}, grad_::AbstractVector{T})
        # ΔA = NaN
        # return ΔA
    end
    return constrain_step, helmholtz_diff_grad!, helmholtz_diff!
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

    # find initial vector for optimizer

    # run optimizer

    # sort phases into foo_1 for gas, foo_2 for liquid

    # return result

    return nothing
end
