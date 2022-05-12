"""
    pressure(substance, nmol, V, RT)

Compute pressure (Pa) of `substance` at given number of moles `nmol` (mol),
total volume `V` (m³) and thermal energy `RT` (J mol⁻¹).
"""
function pressure(substance::AbstractEoSComponent, nmol::Real, V::Real, RT::Real)
    error("NotImplemented")
end

"""
    pressure(mixture, nmol, volume, RT[; buf])

Return pressure (Pa) of `mixture` at given

- `nmol::AbstractVector`: composition (molar parts) or number of moles (mol)
- `volume::Real`: molar volume (m³ mol⁻¹) or volume (m³)
- `RT::Real`: thermal energy (J mol⁻¹)

Allocations may be avoided by passing `buf`.

# Keywords
- `buf::AbstractEoSThermoBuffer`: buffer for intermediate calculations.

See also: [`thermo_buffer`](@ref).
"""
function pressure(
    mixture::AbstractEoSMixture,
    nmol::AbstractVector,
    volume::Real,
    RT::Real;
    buf::AbstractEoSThermoBuffer=thermo_buffer(mixture),
)
    error("NotImplemented")
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

# Seems not necessary for each EoS
function wilson_saturation_pressure(substance::AbstractEoSComponent, RT::Real) end

"""
    compressibility(mixture, nmol, pressure, RT[, phase='g'][; buf])

Compute compressibility (z-factor) of `mixture` in `phase` at given

- composition in moles `nmol` (mol)
- pressure `pressure` (Pa)
- thermal energy `RT` (J mol⁻¹)

# Optional arguments

- `phase::AbstractChar='g'`: specifies phase of `mixture` (`'g'` for gas, `'l'` for liquid)

# Keywords:
- `buf::Union{BrusilovskyThermoBuffer,NamedTuple,AbstractDict}`: buffer for intermediate
    calculations (see [`pressure`](@ref))

See also: [`pressure`](@ref)
"""
function compressibility(
    mix::AbstractEoSMixture,
    nmol::AbstractVector,
    pressure,
    RT,
    phase::AbstractChar;
    buf::AbstractEoSThermoBuffer=thermo_buffer(mix),
)
    zroots = zfactors(mix, nmol, pressure, RT; buf=buf)
    return zfactorchoose(zroots, phase)
end

"Compressibility factors of a `mixture`."
function zfactors(
    mixture::AbstractEoSMixture,
    nmol::AbstractVector,
    pressure,
    RT::Real;
    buf::AbstractEoSThermoBuffer=thermo_buffer(mix),
)
    error("NotImplemented")
end

# TODO: Symbols instead of chars
function zfactorchoose(factors, phase::AbstractChar)
    if phase == 'g'
        return maximum((z) -> isnan(z) ? -Inf : z, factors)
    elseif phase == 'l'
        return minimum((z) -> isnan(z) ? Inf : z, factors)
    else
        throw(DomainError(repr(phase), "Phase must be 'g' or 'l'"))
    end
end
