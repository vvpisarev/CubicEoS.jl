#=
    IdealIdentityState

Variables that make ideal part of hessian over moles to identity matrix.

state = [ 2 √Nᵢ arcsin(√N'ᵢ / √Nᵢ); √N V'/V ], where N = Σᵢ Nᵢ (total moles).
         |-----------------------| |-------|
                    α'ᵢ               α'V
=#

struct IdealIdentityState{V<:AbstractVector} <: AbstractVTFlashState
    x::V
end

function nmolvol(s::IdealIdentityState, nmolb::V, volumeb::T) where {V, T}
    x = value(s)
    αnmol = @view x[1:end-1]
    αvol = x[end]

    nmol = @. nmolb * (sin(αnmol/(2 * sqrt(nmolb))))^2

    volume = (αvol * volumeb) / sqrt(sum(nmolb))
    return (nmol, volume)
end

function IdealIdentityState(
    concentration::AbstractVector,
    saturation::Real,
    nmolb::AbstractVector,
    volumeb::Real,
)
    vol = saturation * volumeb
    nmol = concentration * vol

    x = Vector{Float64}(undef, length(nmolb) + 1)

    @. x[1:end-1] = 2 * sqrt(nmolb) * asin(sqrt(nmol/nmolb))
    x[end] = sqrt(sum(nmolb)) * saturation

    return IdealIdentityState(x)
end

function gradient!(
    grad::AbstractVector,
    state::IdealIdentityState,
    mix::BrusilovskyEoSMixture,
    nmolb::AbstractVector,
    volumeb::Real,
    RT::Real;
    buf::BrusilovskyThermoBuffer=thermo_buffer(mix),
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
    state::IdealIdentityState,
    mix::BrusilovskyEoSMixture,
    nmolb::AbstractVector,
    volumeb::Real,
    RT::Real;
    buf::BrusilovskyThermoBuffer=thermo_buffer(mix),
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

function __vt_flash_optim_closures(
    state1::IdealIdentityState,
    mix::BrusilovskyEoSMixture,
    nmolb::AbstractVector,
    volumeb::Real,
    RT::Real,
)
    state1x = value(state1)
    buf = thermo_buffer(mix)

    helmb = helmholtz(mix, nmolb, volumeb, RT; buf=buf)

    nmol2 = similar(nmolb, Float64)

    # Upper limits for state variables from definition 0 < x[i] < u[i]
    # Necessary for constrain_step
    xupperfromdef = [π * sqrt.(nmolb); sqrt(sum(nmolb))]

    # covolumes = [nmolb .* components(mix).b; -volumeb]

    function clsr_constrain_step(x::AbstractVector, dir::AbstractVector)
        state1x .= x
        α = Inf

        # Finding upper limit
        ## non-negativness: 0 < xi + α di < 1
        for (i, (xi, di, ui)) in enumerate(zip(state1x, dir, xupperfromdef))
            if iszero(di) && !(0 < xi < ui)
                throw(ConstrainStepZeroDirectionError(i, xi))
            end
            α = di < 0 ? min(α, -xi/di) : min(α, (ui - xi)/di)
        end

        # ## covolumes
        # state1xdotcov = dot(state1x, covolumes)
        # dirdotcov = dot(dir, covolumes)
        # if dirdotcov > 0
        #     α = min(α, -state1xdotcov/dirdotcov)
        # elseif dirdotcov < 0
        #     α = min(α, (sum(covolumes) - state1xdotcov)/dirdotcov)
        # end

        # Finding lower limit
        ## Non-negativness
        αlo = -Inf
        for (xi, di, ui) in zip(state1x, dir, xupperfromdef)
            αlo = !iszero(di) && di > 0 ? max(αlo, -xi/di) : max(αlo, (ui - xi)/di)
        end

        # ## covolumes
        # if dirdotcov < 0
        #     αlo = max(αlo, -state1xdotcov/dirdotcov)
        # elseif dirdotcov > 0
        #     αlo = max(αlo, (sum(covolumes)-state1xdotcov)/dirdotcov)
        # end

        α ≤ αlo && throw(ConstrainStepLowerBoundError(x, dir))

        return α
    end

    function clsr_gradient!(grad::AbstractVector, x::AbstractVector)
        state1x .= x
        grad = gradient!(grad, state1, mix, nmolb, volumeb, RT; buf=thermo_buf)
        return grad
    end

    function clsr_helmdiff!(x::AbstractVector, grad::AbstractVector)
        state1x .= x

        nmol1, vol1 = nmolvol(state1, nmolb, volumeb)
        a1 = helmholtz(mix, nmol1, vol1, RT; buf=buf)

        @. nmol2 = nmolb - nmol1
        vol2 = volumeb - vol1
        a2 = helmholtz(mix, nmol2, vol2, RT; buf=buf)

        Δa = (a1 + a2) - helmb
        grad = gradient!(grad, state1, mix, nmolb, volumeb, RT; buf=buf)

        return Δa, grad
    end

    return clsr_constrain_step, clsr_gradient!, clsr_helmdiff!
end
