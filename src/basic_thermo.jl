#=
Basic thermodynamic properties
=#
using LinearAlgebra

include("solvecubic.jl")

@inline ncomponents(mix::BrusilovskyEoSMixture) = length(mix.components)
@inline Base.length(mix::BrusilovskyEoSMixture) = length(mix.components)

"""
    a_coef(substance, RT)

Returns EoS coefficient ``a(T, Ψ) = a_c ϕ(T, Ψ)`` of `substance` at given `RT`.
```math
ϕ(T, psi) = [ 1 + Ψ (1 - T_r^0.5) ]^2
```
Reference: Brusilovsky2002[section: 5.5.2 (algorithm step 3), eq: (4.34), see p.142 and p.164]
"""
function a_coef(substance::BrusilovskyEoSComponent, RT::Real)
    psi = substance.Psi
    phi = ( 1.0 + psi * (1.0 - sqrt(RT / substance.RTc)) )^2
    return substance.ac * phi
end

"""
    pressure(substance, υ, RT)

Computes pressure (Pa) of `substance` at given molar volume `υ` (m³ mol⁻¹) and thermal energy `RT` (J mol⁻¹).
"""
function pressure(substance::BrusilovskyEoSComponent, υ::Real, RT::Real)
    acoeff = a_coef(substance, RT)
    b = substance.b
    c = substance.c
    d = substance.d
    return RT / (υ - b) - acoeff / ((υ + c) * (υ + d))
end

"""
    pressure(substance, nmol, V, RT)

Computes pressure (Pa) of `substance` at given number of moles `nmol` (mol), total volume `V` (m³) and thermal energy `RT` (J mol⁻¹).
"""
function pressure(substance::BrusilovskyEoSComponent, nmol::Real, V::Real, RT::Real)
    acoeff = a_coef(substance, RT)
    b = substance.b
    c = substance.c
    d = substance.d
    υ = V / nmol
    return RT / (υ - b) - acoeff / ((υ + c) * (υ + d))
end

"""
    pressure(substance; nmol = 1, volume, temperature)

Computes pressure (Pa) of `substance` at given number of moles `nmol` (mol), total volume (m³) and temperature (K).
"""
function pressure(substance::BrusilovskyEoSComponent;
                  nmol::Real = 1,
                  volume::Real,
                  temperature::Real)
    RT = GAS_CONSTANT_SI * temperature
    return pressure(substance, nmol, volume, RT)
end

"""
    wilson_saturation_pressure(substance, RT)

Returns approximate saturation pressure of `substance` at `RT` (J mol⁻¹).

Reference: Brusilovsky2002 [p 272, eq 5.4]
"""
function wilson_saturation_pressure(substance::BrusilovskyEoSComponent, RT::Real)
    return wilson_saturation_pressure(
        substance.Pc, substance.RTc, substance.acentric_factor, RT
    )
end

"""
    wilson_saturation_pressure(Pc::Real, RTc::Real, acentric_factor::Real, RT::Real)

Approximate saturation pressure at `RT` using Wilson correlation.

# Arguments
- `Pc::Real` - critical pressure
- `RTc::Real` - gas constant * critical temperature
- `acentric_factor::Real` - acentric factor
- `RT::Real` - gas constant * temperature

Reference: Brusilovsky2002 [p 272, eq 5.4], Mikyska2010 DOI 10.1002/aic.12387 [Algorithm step 1]
"""
function wilson_saturation_pressure(Pc::Real, RTc::Real, acentric_factor::Real, RT::Real)
    return Pc * exp(5.373 * (1.0 + acentric_factor) * (1.0 - RTc / RT))
end

"""
    eos_parameters(mixture::BrusilovskyEoSMixture, nmol, RT[; buf])

Returns `Am`, `Bm`, `Cm`, `Dm` coefficients and `aij` matrix of EoS for `mixture` with
composition `nmol` and thermal energy `RT`. Allocations may be avoided by passing `buf`.

# Arguments
- `nmol::AbstractVector`: Vector of composition
- `RT::Real`: Thermal energy (J mol⁻¹)

# Keywords
- `buf`: see ?pressure(mixture)
"""
function eos_parameters(
    mixture::BrusilovskyEoSMixture{T},
    nmol::AbstractVector{<:Real},
    RT::Real;
    buf = NamedTuple(),
) where {T}
    return __eos_parameters_impl__(mixture, nmol, RT, buf)
end

# When `buf` is BrusilovskyThermoBuffer
function __eos_parameters_impl__(
    mixture::BrusilovskyEoSMixture{T},
    nmol::AbstractVector{<:Real},
    RT::Real,
    buf::BrusilovskyThermoBuffer{T},
) where {T}
    aij = buf.matr
    ai = buf.vec1
    return __eos_parameters_impl__(mixture, nmol, RT, ai, aij)
end

# When `buf` is mapping-like object
function __eos_parameters_impl__(
    mixture::BrusilovskyEoSMixture{T},
    nmol::AbstractVector{<:Real},
    RT::Real,
    buf::Union{NamedTuple, AbstractDict},
) where {T}
    nc = length(nmol)
    aij = haskey(buf, :aij) ? buf[:aij] : Matrix{T}(undef, nc, nc)
    ai = haskey(buf, :ai) ? buf[:ai] : Vector{T}(undef, nc)
    return __eos_parameters_impl__(mixture, nmol, RT, ai, aij)
end

# Core function
# See Brusilovsky2002 [p 187, eqs 4.159, 4.163-4.165] and [p 189, eq 4.168]
function __eos_parameters_impl__(
    mixture::BrusilovskyEoSMixture{T},
    nmol::AbstractVector{<:Real},
    RT::Real,
    auxv::AbstractVector{<:Real},
    auxm::AbstractMatrix{<:Real},
) where {T}
    ai, aij = auxv, auxm

    map!(comp -> a_coef(comp, RT), ai, mixture.components)

    Bm = Cm = Dm = zero(T)
    @inbounds for i in eachindex(nmol, mixture.components)
        Bm += nmol[i] * mixture.components[i].b
        Cm += nmol[i] * mixture.components[i].c
        Dm += nmol[i] * mixture.components[i].d
    end

    temp = RT / GAS_CONSTANT_SI - 273.16
    eij, gij, hij = mixture.eij, mixture.gij, mixture.hij
    aij .= (one(T) .- (eij .+ temp .* (gij .+ temp .* hij))) .* sqrt.(ai .* ai')
    Am = dot(nmol, aij, nmol)  # Am = nmolᵀ aij nmol
    return Am, Bm, Cm, Dm, aij
end

"""
    pressure(mixture, nmol, volume, RT[, buf])

Returns pressure (Pa) of `mixture` at given

- `nmol::AbstractVector`: composition (molar parts) or number of moles (mol)
- `volume::Real`: molar volume (m³ mol⁻¹) or volume (m³)
- `RT::Real`: thermal energy (J mol⁻¹)

Allocations may be avoided by passing `buf`.

# Keywords
- `buf::Union{BrusilovskyThermoBuffer,NamedTuple,AbstractDict}`: Buffer for intermediate
    calculations. In case of `NamedTuple` and `AbstractDict` `buf` should contain `buf[:ai]`
    `NC = ncomponents(mixture)` vector and `buf[:aij]` NCxNC matrix.
"""
function pressure(
    mixture::BrusilovskyEoSMixture,
    nmol::AbstractVector,
    volume::Real,
    RT::Real;
    buf = NamedTuple(),
)
    Am, Bm, Cm, Dm, _ = eos_parameters(mixture, nmol, RT; buf = buf)
    return sum(nmol) * RT / (volume - Bm) - Am/((volume + Cm) * (volume + Dm))
end

"""
    compressibility(mixture, χ, P, RT[, phase='g'])

Computes compressibility (z-factor) of `mixture` in `phase` at given

- composition in moles `χ` (mol)
- pressure `P` (Pa)
- thermal energy `RT` (J mol⁻¹)

# Optional arguments

- `phase::Char='g'` - specifies phase of `mixture` (`'g'` for gas, `'l'` for liquid)
- `buffer::BrusilovskyEoSMixtureBuffer` - buffer for intermediate calculations
"""
function compressibility(
    mix::BrusilovskyEoSMixture,
    χ::AbstractVector,
    P,
    RT,
    phase::AbstractChar='g';
    buf = NamedTuple()
)

    phase in ('g', 'l') || throw(DomainError(repr(phase), "Phase must be 'g' or 'l'"))

    ntotal = sum(χ)
    am, bm, cm, dm = eos_parameters(mix, χ, RT; buf = buf)
    prt = P / (RT * ntotal)
    am *= prt / (RT * ntotal)
    bm *= prt
    cm *= prt
    dm *= prt
    zroots = solve_cubic(
        1.0,
        cm + dm - bm - 1.0,
        am - bm * cm + cm * dm - bm * dm - dm - cm,
        -bm * cm * dm - cm * dm - am * bm
    )

    if phase == 'g'
        return maximum(z for z in zroots if z > 0.0) # NaN is filtered
    else # phase == 'l'
        return minimum(z for z in zroots if z > 0.0)
    end
end
