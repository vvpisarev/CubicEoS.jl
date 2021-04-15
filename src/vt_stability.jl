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
    nc = ncomponents(mix)
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
    ) where {T}

    is_converged = false

    try
        result = DescentMethods.optimize!(optmethod, Dfunc!, η_,
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
    mix::BrusilovskyEoSMixture,
    nmol::AbstractVector,
    volume,
    RT)

    nc = length(mix)
    jacobian_test = zeros(nc, nc)
    aij = zeros(nc, nc)

    loga_parent = similar(nmol)
    loga_test = similar(nmol)

    aux1, aux2 = similar(nmol), similar(nmol)

    log_activity!(loga_parent, aux1, aij, mix, nmol, volume, RT)
    loga_parent .+= log.(nmol ./ volume)
    p_parent = pressure!(aij, aux1, mix, nmol, volume, RT)

    thresh = -1e-5

    b = [comp.b for comp in components(mix)]
    function stabilitytest!(nmol_test, grad)
        log_activity!(grad, aux1, aij, mix, nmol_test, 1, RT)
        p_test = pressure!(aij, aux1, mix, nmol_test, 1, RT)
        grad .+= log.(nmol_test)
        grad .-= loga_parent
        D = dot(grad, nmol_test) - (p_test - p_parent) / RT
        return D, grad
    end

    function stability_step(nmol_test, d)
        αm = convert(eltype(nmol_test), Inf)
        @inbounds for i in eachindex(nmol_test, d)
            α = -nmol_test[i] / d[i]
            if 0 < α < αm
                αm = α
            end
        end
        α = (1 - nmol_test' * b) / (d' * b)  # Σ(xᵢ+αdᵢ)bᵢ < 1, covolume constrain
        if 0 < α < αm
            αm = α
        end
        return αm
    end

    optmethod = DescentMethods.CholBFGS(nmol)

    comp = mix.components
    p_sat = [wilson_saturation_pressure(comp[i], RT) for i in 1:nc]
    # Initial - gas
    p_init = dot(p_sat, nmol) / sum(nmol)
    nmol_test = nmol .* p_sat ./ p_init

    # Test - gas
    z_gg = compressibility(mix, nmol_test, p_init, RT, 'g')
    ntest_gg = nmol_test .* (p_init / (z_gg * RT * sum(nmol_test)))

    _, jacobian_test = log_activity_wj!(loga_test, jacobian_test, aij, aux1, aux2, mix, ntest_gg, 1, RT)
    for i in 1:nc
        jacobian_test[i,i] += 1 / ntest_gg[i]
    end
    DescentMethods.reset!(optmethod, ntest_gg, jacobian_test)

    D_gg = vt_stability_optim_try!(optmethod, ntest_gg, loga_test, stabilitytest!, stability_step)

    if D_gg < thresh
        return false, optmethod.x
    end

    # Test - liquid
    z_gl = compressibility(mix, nmol_test, p_init, RT, 'l')
    ntest_gl = nmol_test .* (p_init / (z_gl * RT * sum(nmol_test)))

    _, jacobian_test = log_activity_wj!(loga_test, jacobian_test, aij, aux1, aux2, mix, ntest_gl, 1, RT)
    for i in 1:nc
        jacobian_test[i,i] += 1 / ntest_gl[i]
    end
    DescentMethods.reset!(optmethod, ntest_gl, jacobian_test)

    D_gl = vt_stability_optim_try!(optmethod, ntest_gl, loga_test, stabilitytest!, stability_step)

    if D_gl < thresh
        return false, optmethod.x
    end

    # Initial - liquid
    nmol_test = nmol ./ p_sat ./ sum(nmol[i] / p_sat[i] for i in 1:nc)
    p_init = dot(p_sat, nmol_test)

    # Test - gas
    
    z_lg = compressibility(mix, nmol_test, p_init, RT, 'g')
    ntest_lg = nmol_test .* (p_init / (z_lg * RT * sum(nmol_test)))

    _, jacobian_test = log_activity_wj!(loga_test, jacobian_test, aij, aux1, aux2, mix, ntest_lg, 1, RT)
    for i in 1:nc
        jacobian_test[i,i] += 1 / ntest_lg[i]
    end
    DescentMethods.reset!(optmethod, ntest_lg, jacobian_test)

    D_lg = vt_stability_optim_try!(optmethod, ntest_lg, loga_test, stabilitytest!, stability_step)

    if D_lg < thresh
        return false, optmethod.x
    end

    # Test - liquid
    z_ll = compressibility(mix, nmol_test, p_init, RT, 'l')
    ntest_ll = nmol_test .* (p_init / (z_ll * RT * sum(nmol_test)))

    _, jacobian_test = log_activity_wj!(loga_test, jacobian_test, aij, aux1, aux2, mix, ntest_ll, 1, RT)
    for i in 1:nc
        jacobian_test[i,i] += 1 / ntest_ll[i]
    end
    DescentMethods.reset!(optmethod, ntest_ll, jacobian_test)

    D_ll = vt_stability_optim_try!(optmethod, ntest_ll, loga_test, stabilitytest!, stability_step)

    if D_ll < thresh
        return false, optmethod.x
    end

    if isnan(D_gg) && isnan(D_gl) && isnan(D_lg) && isnan(D_ll)  # все 4 попытки провалились
        error("VTStability: all tries have failed")
    end

    return true, optmethod.x
end