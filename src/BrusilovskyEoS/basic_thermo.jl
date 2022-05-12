function CubicEoS.pressure(
    substance::BrusilovskyEoSComponent,
    nmol::Real,
    V::Real,
    RT::Real,
)
    acoeff = a_coef(substance, RT)
    b = substance.b
    c = substance.c
    d = substance.d
    υ = V / nmol
    return RT / (υ - b) - acoeff / ((υ + c) * (υ + d))
end

"""
    a_coef(substance, RT)

Return EoS coefficient ``a(T, Ψ) = a_c ϕ(T, Ψ)`` of `substance` at given `RT`.
```math
ϕ(T, Ψ) = [ 1 + Ψ (1 - T_r^0.5) ]^2
```
Reference: Brusilovsky2002[section: 5.5.2 (algorithm step 3), eq: (4.34), see p.142 and p.164]
"""
function a_coef(substance::BrusilovskyEoSComponent, RT::Real)
    psi = substance.Psi
    phi = ( 1.0 + psi * (1.0 - sqrt(RT / substance.RTc)) )^2
    return substance.ac * phi
end

"Reference: Brusilovsky2002 [p 272, eq 5.4]"
function wilson_saturation_pressure(x::BrusilovskyEoSComponent, RT::Real)
    return CubicEoS.wilson_saturation_pressure(x.Pc, x.RTc, x.acentric_factor, RT)
end

function CubicEoS.pressure(
    mixture::BrusilovskyEoSMixture,
    nmol::AbstractVector,
    volume::Real,
    RT::Real;
    buf::BrusilovskyThermoBuffer=thermo_buffer(mixture),
)
    Am, Bm, Cm, Dm, _ = eos_parameters(mixture, nmol, RT; buf=buf)
    return sum(nmol) * RT / (volume - Bm) - Am/((volume + Cm) * (volume + Dm))
end

function CubicEoS.zfactors(
    mix::BrusilovskyEoSMixture,
    molfrac::AbstractVector,
    pressure,
    RT::Real;
    buf::BrusilovskyThermoBuffer=thermo_buffer(mix),
)
    ntotal = sum(molfrac)
    am, bm, cm, dm = eos_parameters(mix, molfrac, RT; buf=buf)
    prt = pressure / (RT * ntotal)
    am *= prt / (RT * ntotal)
    bm *= prt
    cm *= prt
    dm *= prt
    zroots = solve_cubic(
        1.0,
        cm + dm - bm - 1.0,
        am - bm * cm + cm * dm - bm * dm - dm - cm,
        -bm * cm * dm - cm * dm - am * bm,
    )
    return zroots
end

"Am, Bm, Cm, Dm coefficients and aij matrix."
function eos_parameters(
    mixture::BrusilovskyEoSMixture,
    nmol::AbstractVector,
    RT::Real;
    buf::BrusilovskyThermoBuffer=thermo_buffer(mixture),
)
    aij = buf.matr
    ai = buf.vec1
    return __eos_parameters_impl(mixture, nmol, RT, ai, aij)
end

# Core function
# See Brusilovsky2002 [p 187, eqs 4.159, 4.163-4.165] and [p 189, eq 4.168]
function __eos_parameters_impl(
    mixture::BrusilovskyEoSMixture{T},
    nmol::AbstractVector{<:Real},
    RT::Real,
    auxv::AbstractVector{<:Real},
    auxm::AbstractMatrix{<:Real},
) where {T}
    ai, aij = auxv, auxm

    comp = components(mixture)
    psi = comp.Psi
    ac = comp.ac
    @. ai = ac * (1 + psi * (1 - sqrt(RT / comp.RTc)) )^2

    bi, ci, di = comp.b, comp.c, comp.d
    Bm = dot(nmol, bi)
    Cm = dot(nmol, ci)
    Dm = dot(nmol, di)

    tempC = RT / CubicEoS.GAS_CONSTANT_SI - 273.16
    eij, gij, hij = mixture.eij, mixture.gij, mixture.hij
    aij .= (one(T) .- (eij .+ tempC .* (gij .+ tempC .* hij))) .* sqrt.(ai .* ai')
    Am = dot(nmol, aij, nmol)  # Am = nmolᵀ aij nmol
    return Am, Bm, Cm, Dm, aij
end
