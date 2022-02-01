#=
Flash based on physical variables: moles [mol] and volume [m³] of a phase.
=#

struct PhysicalState{V<:AbstractVector} <: AbstractVTFlashState
    x::V
end

nmolvol(s::PhysicalState, nmolb::AbstractVector, volumeb::Real) = (s.x[1:end-1], s.x[end])

function PhysicalState(
    concentration::AbstractVector,
    saturation::Real,
    nmolb::AbstractVector,
    volumeb::Real
)
    x = similar(nmolb, Float64, length(nmolb) + 1)
    volume1 = saturation * volumeb
    @. x[1:end-1] = concentration * volume1
    x[end] = volume1
    return PhysicalState(x)
end

function gradient!(
    grad::AbstractVector,
    state::PhysicalState,
    mix,
    nmolb,
    volumeb,
    RT;
    buf::BrusilovskyThermoBuffer=thermo_buffer(mix),
)
    nmol1 = @view state.x[1:end-1]
    volume1 = state.x[end]

    grad = nvtgradient!(grad, mix, nmol1, volume1, RT; buf=buf)

    grad2 = similar(grad)
    nmol2 = similar(nmol1)
    @. nmol2 = nmolb - nmol1
    volume2 = volumeb - volume1
    grad2 = nvtgradient!(grad2, mix, nmol2, volume2, RT; buf=buf)

    grad .-= grad2

    return grad
end

function hessian!(
    hess::AbstractMatrix,
    state::PhysicalState,
    mix,
    nmolb,
    volumeb,
    RT;
    buf::BrusilovskyThermoBuffer=thermo_buffer(mix),
)
    nmol1 = @view state.x[1:end-1]
    volume1 = state.x[end]

    hess = nvthessian!(hess, mix, nmol1, volume1, RT; buf=buf)

    hess2 = similar(hess)
    nmol2 = similar(nmol1)
    @. nmol2 = nmolb - nmol1
    volume2 = volumeb - volume1
    hess2 = nvthessian!(hess2, mix, nmol2, volume2, RT; buf=buf)

    hess .+= hess2
    return hess
end

function __vt_flash_optim_closures(
    state1::PhysicalState,
    mix::BrusilovskyEoSMixture,
    nmolb::AbstractVector,
    volumeb::Real,
    RT::Real,
)
    thermo_buf = thermo_buffer(mix)
    state1x = value(state1)

    state2 = PhysicalState(similar(state1x))
    state2x = value(state2)

    stateb = PhysicalState(similar(state1x))
    statebx = value(stateb)
    statebx .= [nmolb; volumeb]
    gradb = similar(state1x)
    gradb = nvtgradient!(gradb, mix, nmolb, volumeb, RT; buf=thermo_buf)

    grad1 = similar(gradb)
    grad2 = similar(gradb)

    covolumes = [(c.b for c in components(mix))..., -1.0]
    statebxdotcov = dot(statebx, covolumes)

    function clsr_constrain_step(x::AbstractVector, dir::AbstractVector)
        state1x .= x
        α = Inf
        # Finding upper limit
        ## Non-negativness: 0 < x[i] + α dir[i] < statexb[i]
        @inbounds for i in eachindex(x)
            if iszero(dir[i]) && !(0 < x[i] < statebx[i])
                error("VTFlash: constrain_step. Zero direction $i, but state[$i] = $(x[i])")
            end
            if dir[i] < 0
                α = min(α, -x[i]/dir[i])
            elseif dir[i] > 0
                α = min(α, (statebx[i] - x[i])/dir[i])
            end
        end

        ## Covolume upper limit
        dirdotcov = dot(dir, covolumes)
        state1xterm = dot(state1x, covolumes)/dirdotcov
        if dirdotcov > 0
            α = min(α, -state1xterm)
        elseif dirdotcov < 0
            α = min(α, statebxdotcov/dirdotcov - state1xterm)
        end

        # Checking lower limit
        αlo = -Inf
        ## Non-negativness
        @inbounds for i in eachindex(x)
            if dir[i] < 0
                αlo = max(αlo, (statebx[i] - x[i])/dir[i])
            elseif dir[i] > 0
                αlo = max(αlo, -x[i]/dir[i])
            end
        end

        ## Covolume lower limit
        if dirdotcov < 0
            αlo = max(αlo, -state1xterm)
        elseif dirdotcov > 0
            αlo = max(αlo, statebxdotcov/dirdotcov - state1xterm)
        end

        # Is max α upper than lower limit?
        α ≤ αlo && error("VTFlash: constrain_step. Lower bound not meet.")

        return α
    end

    function clsr_gradient!(grad::AbstractVector, x::AbstractVector)
        state1x .= x
        grad = gradient!(grad, state1, mix, nmolb, volumeb, RT; buf=thermo_buf)
        return grad
    end

    function clsr_helmdiff!(x::AbstractVector, grad::AbstractVector)
        state1x .= x

        nmol1 = @view state1x[1:end-1]
        volume1 = state1x[end]
        grad1 = nvtgradient!(grad1, mix, nmol1, volume1, RT; buf=thermo_buf)

        @. state2x = statebx - state1x
        nmol2 = @view state2x[1:end-1]
        volume2 = state2x[end]
        grad2 = nvtgradient!(grad2, mix, nmol2, volume2, RT; buf=thermo_buf)

        @. grad = grad1 - grad2

        grad2 .-= gradb  # now it is ∇₂ - ∇base
        a = dot(state1x, grad) + dot(statebx, grad2)

        return a, grad
    end

    return clsr_constrain_step, clsr_gradient!, clsr_helmdiff!
end
