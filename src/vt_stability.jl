#=
VT stability algorithm
=#
struct VTStabilityBuffer{B}
    solver_buffer::B
end

@inline function Base.getproperty(b::VTStabilityBuffer, p::Symbol)
    sol_buf = getfield(b, :solver_buffer)
    return getproperty(sol_buf, p)
end

@inline function Base.getindex(b::VTStabilityBuffer, p::Symbol)
    sol_buf = getfield(b, :solver_buffer)
    return getproperty(sol_buf, p)
end

@inline function Base.propertynames(b::VTStabilityBuffer)
    sol_buf = getfield(b, :solver_buffer)
    return propertynames(sol_buf)
end

function VTStabilityBuffer(mix::BrusilovskyEoSMixture)
    thermo_buf = BrusilovskyThermoBuffer(mix)
    loga_parent = similar(thermo_buf.vec1)
    loga_test = similar(loga_parent)
    p_sat = similar(loga_parent)
    bi = similar(loga_parent)
    nmol_test = similar(loga_parent)
    jacobian = similar(thermo_buf.matr)

    solver_buffer = (
        ;
        loga_parent = loga_parent,
        loga_test = loga_test,
        p_sat = p_sat,
        bi = bi,
        nmol_test = nmol_test,
        jacobian = jacobian,
        thermo_buf = thermo_buf,
    )
    return VTStabilityBuffer(solver_buffer)
end

function VTStabilityBuffer(
    thermo_buf::BrusilovskyThermoBuffer
)
    loga_parent = similar(thermo_buf.vec1)
    loga_test = similar(loga_parent)
    p_sat = similar(loga_parent)
    bi = similar(loga_parent)
    nmol_test = similar(loga_parent)
    jacobian = similar(thermo_buf.matr)

    solver_buffer = (
        ;
        loga_parent = loga_parent,
        loga_test = loga_test,
        p_sat = p_sat,
        bi = bi,
        nmol_test = nmol_test,
        jacobian = jacobian,
        thermo_buf = thermo_buf,
    )
    return VTStabilityBuffer(solver_buffer)
end

"""
    vt_stability_buffer(mix)

Create a buffer for intermediate calculations of one-phase stability of a mixture.

For more info see ?vt_stability(mixture).

See also: [`vt_stability`](@ref), [`thermo_buffer`](@ref)
"""
vt_stability_buffer(x) = VTStabilityBuffer(x)

struct VTStabilityResult{T}
    # converged::Bool  # for future
    isstable::Bool
    energy_density::T
    concentration::Vector{T}

    function VTStabilityResult{T}(isstable::Bool, energy_density, concentration) where {T}
        return new{T}(isstable, energy_density, copy(concentration))
    end
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
        result = Downhill.optimize!(Dfunc!, optmethod, η_;
            gtol=1e-3,
            maxiter=1000,
            constrain_step=maxstep,
            reset=false,
        )
        is_converged = result.converged
    catch e
        @warn e
        is_converged = false
    end
    D_min, _ = Dfunc!(optmethod.x, grad_)
    return D_min
end

function vt_stability(
    mix::BrusilovskyEoSMixture{T},
    nmol::AbstractVector,
    volume,
    RT;
    vtstab_buf::VTStabilityBuffer = VTStabilityBuffer(mix),
) where {T}

    nc = length(mix)

    loga_parent = vtstab_buf.loga_parent
    loga_test = vtstab_buf.loga_test
    jacobian = vtstab_buf.jacobian
    thermo_buf = vtstab_buf.thermo_buf
    bi = vtstab_buf.bi
    p_sat = vtstab_buf.p_sat
    nmol_test = vtstab_buf.nmol_test

    log_c_activity!(loga_parent, mix, nmol, volume, RT; buf = thermo_buf)
    loga_parent .+= log.(nmol ./ volume)
    p_parent = pressure(mix, nmol, volume, RT; buf = thermo_buf)

    thresh = -1e-5

    comp = components(mix)
    map!(subst -> subst.b, bi, comp)
    p_sat .= wilson_saturation_pressure.(comp, RT)
    function stabilitytest!(nmol_test, grad)
        log_c_activity!(grad, mix, nmol_test, 1, RT; buf = thermo_buf)
        p_test = pressure(mix, nmol_test, 1, RT; buf = thermo_buf)
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
        _, jacobian = log_c_activity_wj!(
            loga_test,
            jacobian,
            mix,
            nmol_test,
            1,
            RT;
            buf = thermo_buf,
        )
        for i in 1:nc
            jacobian[i,i] += 1 / nmol_test[i]
        end
        Downhill.reset!(optmethod, nmol_test, jacobian)
        return vt_stability_optim_try!(
            optmethod,
            nmol_test,
            loga_test,
            stabilitytest!,
            stability_step,
        )
    end

    optmethod = Downhill.CholBFGS(nmol)
    results = Vector{VTStabilityResult{T}}(undef, 4)

    # Initial - gas
    p_init = dot(p_sat, nmol) / sum(nmol)
    nmol_test .= nmol .* p_sat ./ p_init

    # Test - gas
    z_gg = compressibility(mix, nmol_test, p_init, RT, 'g')
    nmol_test .*= p_init / (z_gg * RT * sum(nmol_test))

    D_gg = D_min(nmol_test, optmethod)
    results[1] = VTStabilityResult{T}(D_gg ≥ thresh, D_gg, optmethod.x)

    # Test - liquid
    z_gl = compressibility(mix, nmol_test, p_init, RT, 'l')
    nmol_test .*= z_gg / z_gl

    D_gl = D_min(nmol_test, optmethod)
    results[2] = VTStabilityResult{T}(D_gl ≥ thresh, D_gl, optmethod.x)

    # Initial - liquid
    nmol_test .= nmol ./ p_sat ./ sum(nmol[i] / p_sat[i] for i in 1:nc)
    p_init = dot(p_sat, nmol_test)

    # Test - gas

    z_lg = compressibility(mix, nmol_test, p_init, RT, 'g')
    nmol_test .*= p_init / (z_lg * RT * sum(nmol_test))

    D_lg = D_min(nmol_test, optmethod)
    results[3] = VTStabilityResult{T}(D_lg ≥ thresh, D_lg, optmethod.x)

    # Test - liquid
    z_ll = compressibility(mix, nmol_test, p_init, RT, 'l')
    nmol_test .*= z_lg / z_ll

    D_ll = D_min(nmol_test, optmethod)
    results[4] = VTStabilityResult{T}(D_ll ≥ thresh, D_ll, optmethod.x)

    if isnan(D_gg) && isnan(D_gl) && isnan(D_lg) && isnan(D_ll)  # все 4 попытки провалились
        error("VTStability: all tries have failed")
    end

    isstable = results[1].isstable && results[2].isstable && results[3].isstable && results[4].isstable

    return isstable, results
end
