include("types.jl")
include("nvt.jl")
include("state_abstract.jl")
include("state_physical.jl")
include("state_idealidentity.jl")

function vt_stability(
    mixture,
    nmol,
    volume,
    RT,
    ::Type{StateVariables};
    buf::AbstractEoSThermoBuffer=thermo_buffer(mixture),
) where {StateVariables<:AbstractVTStabilityState}
    # form base state
    basestate = VTStabilityBaseState(mixture, nmol, volume, RT; buf=buf)

    # prepare TPD closures: TPD, gradient, constrain_step
    tpd_fdf!, tpd_df! = __vt_stability_tpd_closures(StateVariables, basestate; buf=buf)
    constrain_step = __vt_stability_step_closure(StateVariables, basestate.mixture.components.b)

    # TODO: prepare stop criterion closure

    # prepare initial guesses
    trial_concentrations = vt_stability_initials_satpressure(mixture, nmol, RT; buf=buf)

    # run optimizer for each guess
    optmethod = Downhill.CholBFGS(basestate.logconcentration)
    results = map(trial_concentrations) do concentration
        trialstate = fromconcentration(StateVariables, concentration)
        return vt_stability!(trialstate, basestate, optmethod;
            tpd_fdf! = tpd_fdf!,
            constrain_step=constrain_step,
            tpd_thresh=-1e-5,
            maxiter=200,
            buf=buf,
        )
    end

    issuccess = any(x -> x.issuccess, results)
    !issuccess && error("VTStability: all tries have failed")

    isstable = all(x -> x.isstable, results)

    return issuccess, isstable, results
end

function vt_stability!(
    trialstate::AbstractVTStabilityState,
    basestate::VTStabilityBaseState,
    optmethod;
    tpd_fdf!::Function,
    constrain_step::Function,
    tpd_thresh::Real=-1e-5,
    maxiter::Integer=200,
    buf::AbstractEoSThermoBuffer=thermo_buffer(basestate.mixture),
)
    testhessian = let n = ncomponents(basestate.mixture)
        Matrix{Float64}(undef, (n, n))
    end
    testhessian = helmholtztpdhessian!(testhessian, trialstate, basestate; buf=buf)

    teststatex = value(trialstate)
    Downhill.reset!(optmethod, teststatex, testhessian)
    optimresult = Downhill.optimize!(tpd_fdf!, optmethod, teststatex,
        gtol=1e-3,
        # TODO: convcond=convcond,
        maxiter=maxiter,
        constrain_step=constrain_step,
        reset=false,
    )
    # Manually update trialstate
    teststatex .= optimresult.argument

    # Result
    tpd_val = Downhill.fnval(optmethod)
    isstable = tpd_val ≥ -abs(tpd_thresh)
    testconc = concentration(trialstate)
    optim = OptimStats(optimresult.converged, optimresult.iterations, optimresult.calls)

    return __dev_VTStabilityResult(true, isstable, tpd_val, testconc, trialstate, optim)
end

struct OptimStats
    converged::Bool
    iters::Int
    calls::Int
end

struct __dev_VTStabilityResult{T<:Real,S<:AbstractVTStabilityState}
    issuccess::Bool
    isstable::Bool
    energy_density::T
    concentration::Vector{T}
    state::S
    optim::OptimStats
end

function __vt_stability_tpd_closures(
    StateVariables::Type{<:AbstractVTStabilityState},
    basestate::VTStabilityBaseState;
    buf::AbstractEoSThermoBuffer=thermo_buffer(basestate.mixture),
)
    # for internal usage
    trialstate = StateVariables(similar(basestate.logconcentration))
    trialx = value(trialstate)

    function clsr_tpd_fdf!(x::AbstractVector, g::AbstractVector)
        trialx .= x
        tpd, g = helmholtztpdwgradient!(g, trialstate, basestate; buf=buf)
        return tpd, g
    end

    function clsr_tpd_df!(g::AbstractVector, x::AbstractVector)
        trialx .= x
        g = helmholtztpdgradient!(g, trialstate, basestate; buf=buf)
        return g
    end

    return clsr_tpd_fdf!, clsr_tpd_df!
end

function __vt_stability_step_closure(
    ::Type{StateVariables},
    covolumes::AbstractVector,
) where {StateVariables<:AbstractVTStabilityState}
    function clsr(x::AbstractVector, direction::AbstractVector)
        return __constrain_step(StateVariables, x, direction, covolumes)
    end
    return clsr
end

# TODO: Check if zfactors are distinct
# TODO: DRY
function vt_stability_initials_satpressure(
    mixture::AbstractEoSMixture,
    nmolbase::AbstractVector,
    RT::Real;
    buf::AbstractEoSThermoBuffer=thermo_buffer(mixture),
)
    psat = wilson_saturation_pressure.(mixture.components, RT)

    # Base is liquid, trial is vapor/liquid
    p_baseliquid = dot(psat, nmolbase) / sum(nmolbase)
    molfrac_trialvapor = (psat ./ p_baseliquid) .* (nmolbase ./ sum(nmolbase))
    zroots_baseliquid = zfactors(mixture, molfrac_trialvapor, p_baseliquid, RT; buf=buf)

    conc_trial_lv = let z = zfactorchoose(zroots_baseliquid, 'g')
        p_baseliquid .* molfrac_trialvapor ./ (z * RT)
    end
    conc_trial_ll = let z = zfactorchoose(zroots_baseliquid, 'l')
        p_baseliquid .* molfrac_trialvapor ./ (z * RT)
    end

    # Base is vapor, trial is vapor/liquid
    molfrac_trialliquid = (nmolbase ./ psat) ./ sum(nmolbase ./ psat)
    p_basevapor = dot(psat, molfrac_trialliquid)
    zroots_basevapor = zfactors(mixture, molfrac_trialliquid, p_basevapor, RT; buf=buf)

    conc_trial_vv = let z = zfactorchoose(zroots_basevapor, 'g')
        p_basevapor .* molfrac_trialliquid ./ (z * RT)
    end
    conc_trial_vl = let z = zfactorchoose(zroots_basevapor, 'l')
        p_basevapor .* molfrac_trialliquid ./ (z * RT)
    end

    return (conc_trial_lv, conc_trial_ll, conc_trial_vv, conc_trial_vl)
end

#= OLD code goes down =#

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
    issuccess::Bool
    isstable::Bool
    energy_density::T
    concentration::Vector{T}

    function VTStabilityResult{T}(issuccess::Bool, isstable::Bool, energy_density, concentration) where {T}
        return new{T}(issuccess, isstable, energy_density, copy(concentration))
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
    isfinished = false
    try
        optresult = Downhill.optimize!(Dfunc!, optmethod, η_;
            gtol=1e-3,
            maxiter=1000,
            constrain_step=maxstep,
            reset=false,
        )
        isfinished = true
    catch e
        if e isa InterruptException
            rethrow(e)
        else
            @warn "$(sprint(showerror, e))"
        end
    end

    Dtry = isfinished ? Downhill.fnval(optmethod) : NaN

    return Dtry
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
    results[1] = VTStabilityResult{T}(!isnan(D_gg), D_gg ≥ thresh, D_gg, optmethod.x)

    # Test - liquid
    z_gl = compressibility(mix, nmol_test, p_init, RT, 'l')
    nmol_test .*= z_gg / z_gl

    D_gl = D_min(nmol_test, optmethod)
    results[2] = VTStabilityResult{T}(!isnan(D_gl), D_gl ≥ thresh, D_gl, optmethod.x)

    # Initial - liquid
    nmol_test .= nmol ./ p_sat ./ sum(nmol[i] / p_sat[i] for i in 1:nc)
    p_init = dot(p_sat, nmol_test)

    # Test - gas

    z_lg = compressibility(mix, nmol_test, p_init, RT, 'g')
    nmol_test .*= p_init / (z_lg * RT * sum(nmol_test))

    D_lg = D_min(nmol_test, optmethod)
    results[3] = VTStabilityResult{T}(!isnan(D_gl), D_lg ≥ thresh, D_lg, optmethod.x)

    # Test - liquid
    z_ll = compressibility(mix, nmol_test, p_init, RT, 'l')
    nmol_test .*= z_lg / z_ll

    D_ll = D_min(nmol_test, optmethod)
    results[4] = VTStabilityResult{T}(!isnan(D_gl), D_ll ≥ thresh, D_ll, optmethod.x)

    issuccess = any(x -> x.issuccess, results)
    !issuccess && error("VTStability: all tries have failed")

    isstable = !any(x -> x.issuccess && !x.isstable, results)

    return issuccess, isstable, results
end
