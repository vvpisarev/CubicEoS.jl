#=
Flash based on fractional (ratio) variables

        [ N'1       N'n   V' ]
state = [ ---, ..., ---, --- ]
        [  N1        Nn   V  ]
=#

struct RatioState{V<:AbstractVector} <: AbstractVTFlashState
    x::V
end

function nmolvol(s::RatioState, nmolb::V, volumeb::T) where {V, T}
    x = value(s)
    nmol1 = nmolb .* x[1:end-1]
    volume1 = volumeb * x[end]
    return (nmol1, volume1)
end

function RatioState(
    concentration::AbstractVector,
    saturation::Real,
    nmolb::AbstractVector,
    volumeb::Real
)
    x = similar(nmolb, Float64, length(nmolb) + 1)
    @. x[1:end-1] = volumeb * saturation * concentration / nmolb
    x[end] = saturation
    return RatioState(x)
end

function gradient!(
    grad::AbstractVector,
    state::RatioState,
    mix::BrusilovskyEoSMixture,
    nmolb::AbstractVector,
    volumeb::Real,
    RT::Real;
    buf::BrusilovskyThermoBuffer=thermo_buffer(mix),
)
    nmol1, vol1 = nmolvol(state, nmolb, volumeb)
    grad = nvtgradient!(grad, mix, nmol1, vol1, RT; buf=buf)

    nmol2 = nmolb - nmol1  # alloc
    vol2 = volumeb - vol1
    grad2 = similar(grad)  # alloc
    grad2 = nvtgradient!(grad2, mix, nmol2, vol2, RT; buf=buf)

    grad .-= grad2
    grad[1:end-1] .*= nmolb
    grad[end] *= volumeb

    return grad
end

function hessian!(
    hess::AbstractMatrix,
    state::RatioState,
    mix::BrusilovskyEoSMixture,
    nmolb::AbstractVector,
    volumeb::Real,
    RT::Real;
    buf::BrusilovskyThermoBuffer=thermo_buffer(mix),
)
    nmol1, vol1 = nmolvol(state, nmolb, volumeb)
    hess = nvthessian!(hess, mix, nmol1, vol1, RT; buf=buf)

    nmol2 = nmolb - nmol1  # alloc
    vol2 = volumeb - vol1
    hess2 = similar(hess)  # alloc
    hess2 = nvthessian!(hess2, mix, nmol2, vol2, RT; buf=buf)

    hess .+= hess2

    # Also, kron and gemm! is applicable
    scale = [nmolb; volumeb] * [nmolb; volumeb]'  # alloc
    hess .*= scale

    return hess
end

function __vt_flash_optim_closures(
    state1::RatioState,
    mix::BrusilovskyEoSMixture,
    nmolb::AbstractVector,
    volumeb::Real,
    RT::Real,
)
    state1x = value(state1)
    buf = thermo_buffer(mix)

    helmb = helmholtz(mix, nmolb, volumeb, RT; buf=buf)

    nmol2 = similar(nmolb, Float64)

    covolumes = [nmolb .* components(mix).b; -volumeb]

    function clsr_constrain_step(x::AbstractVector, dir::AbstractVector)
        state1x .= x
        α = Inf

        # Finding upper limit
        ## non-negativness: 0 < xi + α di < 1
        for (i, (xi, di)) in enumerate(zip(state1x, dir))
            if iszero(di) && !(0 < xi < 1)
                throw(ConstrainStepZeroDirectionError(i, xi))
            end
            α = di < 0 ? min(α, -xi/di) : min(α, (1 - xi)/di)
        end

        ## covolumes
        state1xdotcov = dot(state1x, covolumes)
        dirdotcov = dot(dir, covolumes)
        if dirdotcov > 0
            α = min(α, -state1xdotcov/dirdotcov)
        elseif dirdotcov < 0
            α = min(α, (sum(covolumes) - state1xdotcov)/dirdotcov)
        end

        # Finding lower limit
        ## Non-negativness
        αlo = -Inf
        for (xi, di) in zip(state1x, dir)
            αlo = !iszero(di) && di > 0 ? max(αlo, -xi/di) : max(αlo, (1-xi)/di)
        end

        ## covolumes
        if dirdotcov < 0
            αlo = max(αlo, -state1xdotcov/dirdotcov)
        elseif dirdotcov > 0
            αlo = max(αlo, (sum(covolumes)-state1xdotcov)/dirdotcov)
        end

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
