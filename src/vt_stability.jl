#=
VT stability algorithm
=#
export vt_stability

using DescentMethods

struct StabilityBuffer{M, P}
    mixture_eos::M
    parent_phase::P
    test_phase::P
end

function StabilityBuffer(mix::BrusilovskyEoSMixture{T}) where {T}
    parent_phase = ()
    test_phase = ()
    return StabilityBuffer(mix, parent_phase, test_phase)
end

"""
Оптимизационная процедура для алгоритма проверки стабильности.
Минимизирует D(η').

Возвращает минимум функции D.
Обновляет аргумент функции D - концентрацию η_.

-> η_, D_thresh, Dfunc!, Dgrad!, isvalid_η, buffer
<- D_min (NaN if not converged)
"""
function vt_stability_optim_try!(
    optmethod,
    η_::AbstractVector,  # вектор концентраций, должен быть инициализирован начальными значениями
    grad_::AbstractVector,
    Dfunc!::Function,  # функция D вида D!(η', ∇D_) -> D::Number
    maxstep::Function,
)

    is_converged = false

    try
        result = DescentMethods.optimize!(
            optmethod,
            Dfunc!,
            η_,
            gtol=1e-3,
            maxiter=1000,
            constrain_step=maxstep,
            reset=false
        )
        η_ .= result.argument

        # is_converged = true
        is_converged = result.converged
    catch e
        @warn e
        is_converged = false
    end

    #if is_converged
        D_min, _ = Dfunc!(η_, grad_) # current_grad взят, чтобы не аллоцировать
    #else
    #    D_min = NaN
    #end
    return D_min
end

function vt_stability(
    mix::BrusilovskyEoSMixture{T},
    nmol::AbstractVector,
    volume,
    RT;
    kwargs...
) where {T}

    nc = length(mix)
    buffers = kwargs.data

    loga_parent = haskey(buffers, :loga_parent) ? buffers[:loga_parent] : Vector{T}(undef, nc)
    loga_test = haskey(buffers, :loga_test) ? buffers[:loga_test] : Vector{T}(undef, nc)
    jacobian = haskey(buffers, :jacobian) ? buffers[:jacobian] : Matrix{T}(undef, nc, nc)
    aux1 = haskey(buffers, :aux1) ? buffers[:aux1] : Vector{T}(undef, nc)
    aux2 = haskey(buffers, :aux2) ? buffers[:aux2] : Vector{T}(undef, nc)
    bi = haskey(buffers, :bi) ? buffers[:bi] : Vector{T}(undef, nc)
    p_sat = haskey(buffers, :p_sat) ? buffers[:p_sat] : Vector{T}(undef, nc)
    nmol_test = haskey(buffers, :nmol_test) ? buffers[:nmol_test] : Vector{T}(undef, nc)

    aij = jacobian
    log_activity(mix, nmol, volume, RT; log_a = loga_parent, ai = aux1, aij = aij)
    loga_parent .+= log.(nmol ./ volume)
    p_parent = pressure(mix, nmol, volume, RT; aij = aij, ai = aux1)

    thresh = -1e-5

    comp = components(mix)
    map!(subst -> subst.b, bi, comp)
    map!(subst -> wilson_saturation_pressure(subst, RT), p_sat, comp)
    function stabilitytest!(nmol_test, grad)
        log_activity(mix, nmol_test, 1, RT; log_a = grad, ai = aux1, aij = aij)
        p_test = pressure(mix, nmol_test, 1, RT; aij = aij, ai = aux1)
        grad .+= log.(nmol_test)
        grad .-= loga_parent
        D = dot(grad, nmol_test) - (p_test - p_parent) / RT
        return D, grad
    end

    function stability_step(nmol_test, dir)
        αm = convert(eltype(nmol_test), Inf)
        @inbounds for i in eachindex(nmol_test, dir)
            α = -nmol_test[i] / dir[i]
            if 0 < α < αm
                αm = α
            end
        end
        α = (1 - nmol_test' * bi) / (dir' * bi)  # Σ(xᵢ+αdᵢ)bᵢ < 1, covolume constrain
        if 0 < α < αm
            αm = α
        end
        return αm
    end

    function D_min(nmol_test, optmethod)
        _, jacobian = log_activity_wj(
            mix,
            nmol_test,
            1,
            RT;
            log_a = loga_test,
            jacobian = jacobian,
            aux1 = aux1,
            aux2 = aux2
        )
        for i in 1:nc
            jacobian[i,i] += 1 / nmol_test[i]
        end
        DescentMethods.reset!(optmethod, nmol_test, jacobian)
        return vt_stability_optim_try!(
            optmethod,
            nmol_test,
            loga_test,
            stabilitytest!,
            stability_step
        )
    end

    optmethod = DescentMethods.CholBFGS(nmol)

    # Initial - gas
    p_init = dot(p_sat, nmol) / sum(nmol)
    nmol_test .= nmol .* p_sat ./ p_init

    # Test - gas
    z_gg = compressibility(mix, nmol_test, p_init, RT, 'g')
    nmol_test .*= p_init / (z_gg * RT * sum(nmol_test))

    D_gg = D_min(nmol_test, optmethod)

    if D_gg < thresh
        return false, optmethod.x
    end

    # Test - liquid
    z_gl = compressibility(mix, nmol_test, p_init, RT, 'l')
    nmol_test .*= z_gg / z_gl

    D_gl = D_min(nmol_test, optmethod)

    if D_gl < thresh
        return false, optmethod.x
    end

    # Initial - liquid
    nmol_test .= nmol ./ p_sat ./ sum(nmol[i] / p_sat[i] for i in 1:nc)
    p_init = dot(p_sat, nmol_test)

    # Test - gas

    z_lg = compressibility(mix, nmol_test, p_init, RT, 'g')
    nmol_test .*= p_init / (z_lg * RT * sum(nmol_test))

    D_lg = D_min(nmol_test, optmethod)

    if D_lg < thresh
        return false, optmethod.x
    end

    # Test - liquid
    z_ll = compressibility(mix, nmol_test, p_init, RT, 'l')
    nmol_test .*= z_lg / z_ll

    D_ll = D_min(nmol_test, optmethod)

    if D_ll < thresh
        return false, optmethod.x
    end

    if isnan(D_gg) && isnan(D_gl) && isnan(D_lg) && isnan(D_ll)  # все 4 попытки провалились
        error("VTStability: all tries have failed")
    end

    return true, optmethod.x
end
