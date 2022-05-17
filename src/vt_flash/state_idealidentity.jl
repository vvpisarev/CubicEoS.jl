"""
    VTFlashIdealIdentityState(x)
    VTFlashIdealIdentityState(concentration, saturation, nmolbase, volumebase)

Variables that make ideal part of Hessian over moles to identity matrix.

**Definition**

```
state = [ 2 √Nᵢ arcsin(√N'ᵢ / √Nᵢ); √N V'/V ]
         |-----------------------| |-------|
                    α'ᵢ               α'V
```
where superscript `'` relates to a trial phase and `N = Σᵢ Nᵢ` is total moles.
"""
struct VTFlashIdealIdentityState{V<:AbstractVector} <: AbstractVTFlashState
    x::V
end

function nmolvol!(nmol, s::VTFlashIdealIdentityState, nmolbase, volumebase)
    x = value(s)
    αnmol = @view x[1:end-1]
    αvol = x[end]

    @. nmol = nmolbase * (sin(αnmol/(2 * sqrt(nmolbase))))^2

    volume = (αvol * volumebase) / sqrt(sum(nmolbase))
    return nmol, volume
end

function VTFlashIdealIdentityState{V}(
    concentration::AbstractVector,
    saturation::Real,
    nmolb::AbstractVector,
    volumeb::Real,
) where {V}
    vol = saturation * volumeb
    nmol = concentration * vol

    x = Vector{Float64}(undef, length(nmolb) + 1)

    @. x[1:end-1] = 2 * sqrt(nmolb) * asin(sqrt(nmol/nmolb))
    x[end] = sqrt(sum(nmolb)) * saturation

    return VTFlashIdealIdentityState{V}(x)
end

@inline VTFlashIdealIdentityState(c, s, n, v) = VTFlashIdealIdentityState{Vector{Float64}}(c, s, n, v)

function gradient!(
    grad::AbstractVector,
    state::VTFlashIdealIdentityState,
    mix::AbstractEoSMixture,
    nmolb::AbstractVector,
    volumeb::Real,
    RT::Real;
    buf::AbstractEoSThermoBuffer=thermo_buffer(mix),
)
    nmol1, vol1 = nmolvol(state, nmolb, volumeb)
    grad = nvtgradient!(grad, mix, nmol1, vol1, RT; buf=buf)

    nmol2 = nmolb - nmol1
    vol2 = volumeb - vol1
    grad2 = similar(grad)
    grad2 = nvtgradient!(grad2, mix, nmol2, vol2, RT; buf=buf)

    grad .-= grad2

    # Scaling
    @. grad[1:end-1] *= sqrt(nmol1 * nmol2 / nmolb)
    grad[end] *= volumeb / sqrt(sum(nmolb))

    return grad
end

function hessian!(
    hess::AbstractMatrix,
    state::VTFlashIdealIdentityState,
    mix::AbstractEoSMixture,
    nmolb::AbstractVector,
    volumeb::Real,
    RT::Real;
    buf::AbstractEoSThermoBuffer=thermo_buffer(mix),
)
    # Here is version that scales nvt-hessian
    # but a direct calculation may be considerd

    nmol1, vol1 = nmolvol(state, nmolb, volumeb)
    grad1 = similar(value(state))
    grad1, hess = nvtgradienthessian!(grad1, hess, mix, nmol1, vol1, RT; buf=buf)

    nmol2, vol2 = nmolb - nmol1, volumeb - vol1
    grad2 = similar(value(state))
    hess2 = similar(hess)
    grad2, hess2 = nvtgradienthessian!(grad2, hess2, mix, nmol2, vol2, RT; buf=buf)

    hess .+= hess2  # now, `hess` = hessian of a' + hessian of a''

    # dNdN part scaling
    @inbounds for j in 1:size(hess, 2)-1, i in 1:size(hess, 1)-1
        hess[i, j] *= sqrt(nmol1[i] * nmol2[i] * nmol1[j] * nmol2[j] / (nmolb[i] * nmolb[j]))
        if i == j
            hess[i, j] += (grad1[i] - grad2[i]) * (nmol2[i] - nmol1[i]) / (2*nmolb[i])
        end
    end

    # dNdV part scaling
    @inbounds for i in 1:size(hess, 1)-1
        sumnmolb = sum(nmolb)
        hess[i, end] *= volumeb * sqrt(nmol1[i]*nmol2[i] / (nmolb[i] * sumnmolb))
    end
    hess[end, 1:end-1] .= hess[1:end-1, end]

    # dV² part
    hess[end, end] *= volumeb^2 / sum(nmolb)

    return hess
end

function physical_constrain_step_closure(
    ::Type{<:VTFlashIdealIdentityState},
    nmolbase::AbstractVector,
    volumebase::Real=NaN,  # ignored
)
    # Upper limits for variables
    xuplims = [π * sqrt.(nmolbase); sqrt(sum(nmolbase))]

    function clsr(
        ::Type{<:VTFlashIdealIdentityState},
        x::AbstractVector,
        direction::AbstractVector,
    )
        αmax = Inf
        for (i, (xi, di, ui)) in enumerate(zip(x, direction, xuplims))
            if iszero(di) && !(0 < xi < ui)
                throw(ConstrainStepZeroDirectionError(i, xi))
            end
            αmax = di < 0 ? min(αmax, -xi/di) : min(αmax, (ui - xi)/di)
        end

        αlow = -Inf
        for (xi, di, ui) in zip(x, direction, xuplims)
            # Zero direction is checked above
            αlow = di > 0 ? max(αlow, -xi/di) : max(αlow, (ui - xi)/di)
        end
        return αlow, αmax
    end
    return clsr
end

# TODO: deprecate, but it stores linear approximation for cubic eos constraint
# function __vt_flash_optim_closures(
#     state1::VTFlashIdealIdentityState,
#     mix::AbstractEoSMixture,
#     nmolb::AbstractVector,
#     volumeb::Real,
#     RT::Real,
# )
#     state1x = value(state1)
#     buf = thermo_buffer(mix)

#     helmb = helmholtz(mix, nmolb, volumeb, RT; buf=buf)

#     nmol2 = similar(nmolb, Float64)

#     # Upper limits for state variables from definition 0 < x[i] < u[i]
#     # Necessary for constrain_step
#     xupperfromdef = [π * sqrt.(nmolb); sqrt(sum(nmolb))]

#     BASECOVDIFF = dot(nmolb, components(mix).b) - volumeb

#     "Check of covolume constrain for phases."
#     function iscovmeet(x::AbstractVector, α::Real, d::AbstractVector)
#         phaseterm = zero(Float64)
#         for (xᵢ, dᵢ, nᵢ, bᵢ) in zip(x[1:end-1], d[1:end-1], nmolb, components(mix).b)
#             phaseterm += nᵢ * bᵢ * sin((xᵢ + α*dᵢ)/(2*sqrt(nᵢ)))^2
#         end
#         volterm = volumeb * (x[end] + α * d[end]) / sqrt(sum(nmolb))
#         cond1 = phaseterm ≤ volterm
#         cond2 = BASECOVDIFF - volumeb ≤ phaseterm - volterm
#         return cond1 && cond2
#     end

#     "Linear approximation to covolume condition for phase 1: `constant + α linear ≤ 0`."
#     function covcondfactors(x::AbstractVector, d::AbstractVector)
#         constant = zero(Float64)
#         linear = zero(Float64)

#         for (xᵢ, dᵢ, nᵢ, bᵢ) in zip(x[1:end-1], d[1:end-1], nmolb, components(mix).b)
#             sqrtnᵢ = sqrt(nᵢ)
#             constant += nᵢ * bᵢ * sin(xᵢ/(2*sqrtnᵢ))^2
#             linear += sqrtnᵢ * bᵢ * sin(xᵢ/sqrtnᵢ) * dᵢ
#         end
#         constant -= volumeb / sqrt(sum(nmolb)) * x[end]
#         linear /= 2
#         linear -= volumeb * d[end] / sqrt(sum(nmolb))

#         return constant, linear
#     end

#     function clsr_constrain_step(x::AbstractVector, dir::AbstractVector)
#         state1x .= x
#         α = Inf

#         # Finding upper limit
#         ## non-negativness: 0 < xi + α di < 1
#         for (i, (xi, di, ui)) in enumerate(zip(state1x, dir, xupperfromdef))
#             if iszero(di) && !(0 < xi < ui)
#                 throw(ConstrainStepZeroDirectionError(i, xi))
#             end
#             α = di < 0 ? min(α, -xi/di) : min(α, (ui - xi)/di)
#         end

#         ## covolumes
#         covconst, covlin = covcondfactors(x, dir)
#         α = covlin > 0 ? min(α, -covconst/covlin) : min(α, (BASECOVDIFF - covconst)/covlin)

#         # Finding lower limit
#         ## Non-negativness
#         αlo = -Inf
#         for (xi, di, ui) in zip(state1x, dir, xupperfromdef)
#             αlo = !iszero(di) && di > 0 ? max(αlo, -xi/di) : max(αlo, (ui - xi)/di)
#         end
#         αlo = max(zero(αlo), αlo)
#         ## covolumes
#         αlo = covlin < 0 ? max(αlo, -covconst/covlin) : max(αlo, (BASECOVDIFF - covconst)/covlin)

#         α ≤ αlo && throw(ConstrainStepLowerBoundError(x, dir))

#         # Check that α suitable for original covolume constrain
#         # Actually, sometimes do not work, but optimization goes OK.
#         # Perhaps, because of bad approximation of covolume constrain,
#         # and α is limitted by non-negativness, which gives absence of second phase.
#         # !iscovmeet(x, α, dir) && error("ConstrainStep: step not found")

#         return α
#     end

#     function clsr_gradient!(grad::AbstractVector, x::AbstractVector)
#         state1x .= x
#         grad = gradient!(grad, state1, mix, nmolb, volumeb, RT; buf=buf)
#         return grad
#     end

#     function clsr_helmdiff!(x::AbstractVector, grad::AbstractVector)
#         state1x .= x

#         nmol1, vol1 = nmolvol(state1, nmolb, volumeb)
#         a1 = helmholtz(mix, nmol1, vol1, RT; buf=buf)

#         @. nmol2 = nmolb - nmol1
#         vol2 = volumeb - vol1
#         a2 = helmholtz(mix, nmol2, vol2, RT; buf=buf)

#         Δa = (a1 + a2) - helmb
#         grad = gradient!(grad, state1, mix, nmolb, volumeb, RT; buf=buf)

#         return Δa, grad
#     end

#     return clsr_constrain_step, clsr_gradient!, clsr_helmdiff!
# end
