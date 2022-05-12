"""
    log_c_activity(mixture, nmol, volume, RT[; buf])

Return vector of ln(c_a) - logarithm of activity coefficients
for components of `mixture` at given `nmol`, `volume`, `RT`.
If buffers are provided as keyword arguments, their contents
are modified during the intermediate calculations.

# Arguments

- `mix`: mixture
- `nmol::AbstractVector`: amount of each component (mol)
- `volume`: volume of the mixture (m³)
- `RT`: thermal energy (J/mol)

# Keywords

- `buf::AbstractThermoBuffer`: buffers for intermediate calculations

# Returns

- `AbstractVector`: the logarithms of activity coefficients of the components at given
    number of moles, volume and temperature

See (Jirí Mikyska, Abbas Firoozabadi // 10.1002/aic.12387)

See also: [`log_c_activity!`](@ref), [`log_c_activity_wj`](@ref), [`log_c_activity_wj!`](@ref)
"""
function log_c_activity(
    mix::AbstractEoSMixture{T},
    nmol::AbstractVector,
    volume::Real,
    RT::Real;
    buf::AbstractEoSThermoBuffer=thermo_buffer(mix),
) where {T}
    nc = ncomponents(mix)
    log_ca = Vector{T}(undef, nc)
    return log_c_activity!(log_ca, mix, nmol, volume, RT; buf=buf)
end

"""
    log_c_activity!(log_ca, mixture, nmol, volume, RT[; buf])

Return vector of ln(c_a) - logarithm of activity coefficients
for components of `mixture` at given `nmol`, `volume`, `RT`.
If buffers are provided as keyword arguments, their contents
are modified during the intermediate calculations. The answer is stored
in `log_ca`.

# Arguments

- `log_ca::AbstractVector`: buffer to store the result
- `mix::BrusilovskyEoSMixture`: mixture
- `nmol::AbstractVector`: amount of each component (mol)
- `volume`: volume of the mixture (m³)
- `RT`: thermal energy (J/mol)

# Keywords

- `buf::AbstractThermoBuffer`: buffers for intermediate calculations

# Returns

- `AbstractVector`: the logarithms of activity coefficients of the components at given
    number of moles, volume and temperature

See (Jirí Mikyska, Abbas Firoozabadi // 10.1002/aic.12387)

See also: [`log_c_activity`](@ref), [`log_c_activity_wj`](@ref), [`log_c_activity_wj!`](@ref)
"""
function log_c_activity!(
    log_ca::AbstractVector,
    mix::AbstractEoSMixture,
    nmol::AbstractVector,
    volume,
    RT;
    buf::AbstractEoSThermoBuffer=thermo_buffer(mix),
)
    error("NotImplemented")
end

"""
    log_c_activity_wj(mixture, nmol, volume, RT[; buf])

Return vector of ln(c_a) - logarithm of activity coefficient
for components of `mixture` at given `nmol`, `volume`, `RT` -
and the jacobian ∂ln(c_a[i]) / ∂n[j].
If buffers are provided as keyword arguments, their contents
are modified during the intermediate calculations.

# Arguments

- `mix`: mixture
- `nmol::AbstractVector`: amount of each component (mol)
- `volume`: volume of the mixture (m³)
- `RT`: thermal energy (J/mol)

# Keywords

- `buf::AbstractThermoBuffer`: buffers for intermediate calculations

# Returns

- `Tuple{AbstractVector,AbstractMatrix}`: the logarithms of activity coefficients
of the components at given number of moles, volume and temperature, and the jacobian
matrix ∂ln(c_a[i]) / ∂n[j].

See also: [`log_c_activity`](@ref), [`log_c_activity!`](@ref), [`log_c_activity_wj!`](@ref)
"""
function log_c_activity_wj(
    mix::AbstractEoSMixture{T},
    nmol::AbstractVector,
    volume::Real,
    RT::Real;
    buf::AbstractEoSThermoBuffer=thermo_buffer(mix),
) where {T}
    nc = ncomponents(mix)
    log_ca = similar(nmol, T)
    jacobian = Matrix{T}(undef, nc, nc)
    return log_c_activity_wj!(log_ca, jacobian, mix, nmol, volume, RT; buf=buf)
end

"""
    log_activity_wj!(log_ca, jacobian, mixture, nmol, volume, RT[; buf])

Return vector of ln(c_a) - logarithm of activity coefficient
for components of `mixture` at given `nmol`, `volume`, `RT` -
and the jacobian ∂ln(c_a[i]) / ∂n[j]. The first two arguments get overwritten
by the result.
If buffers are provided as keyword arguments, their contents
are modified during the intermediate calculations.

# Arguments

- `log_ca::AbstractVector`: a vector to store the activity coefficients
- `jacobian::AbstractMatrix`: a matrix to store ∂ln(c_a[i]) / ∂n[j]
- `mix`: mixture
- `nmol::AbstractVector`: amount of each component (mol)
- `volume`: volume of the mixture (m³)
- `RT`: thermal energy (J/mol)

# Keywords

- `buf::BrusilovskyThermoBuffer`: buffers for intermediate calculations
    (see [`thermo_buffer`](@ref))

# Returns

- `Tuple{AbstractVector,AbstractMatrix}`: the logarithms of activity coefficients
of the components at given number of moles, volume and temperature, and the jacobian
matrix ∂ln(c_a[i]) / ∂n[j]. The values are aliases for the first two arguments.

See also: [`log_c_activity`](@ref), [`log_c_activity!`](@ref), [`log_c_activity_wj`](@ref)
"""
function log_c_activity_wj!(
    log_ca::AbstractVector,
    jacobian::AbstractMatrix,
    mix::AbstractEoSMixture,
    nmol::AbstractVector,
    volume::Real,
    RT::Real;
    buf::AbstractEoSThermoBuffer=thermo_buffer(mix),
)
    error("NotImplemented")
end
