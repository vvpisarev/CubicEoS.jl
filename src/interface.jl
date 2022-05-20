#
# Interface of a component
#

"""
    molar_mass(c::AbstractEoSComponent)

**Required**. Molar mass of a component [kg mole⁻¹].
"""
molar_mass(c::AbstractEoSComponent) = error("NotImplemented")

"""
    carbon_number(c::AbstractEoSComponent)

**Required**. Number of carbons atoms in the hydrocarbon chain.
"""
carbon_number(c::AbstractEoSComponent) = error("NotImplemented")

"""
    name(c::AbstractEoSComponent)
    name(m::AbstractEoSMixture)

**Required**. Name of component (mixture) as `String`-like object.
By default, name of mixture is concatenation of its components names.
"""
name

name(c::AbstractEoSComponent) = error("NotImplemented")

"""
    load(::Type{<:AbstractEoSComponent}; name::AbstractString, databases)

Load component `::T` by its `name` by joining the data on this species in `databases`.
"""
load(::Type{T}; name, component_dbs) where {T<:AbstractEoSComponent} = error("NotImplemented")

#
# Interface of a mixture
#

"""
    ncomponents(mixture::AbstractEoSMixture) -> Int

**Required**. Number of components in `mixture`.
"""
ncomponents(m::AbstractEoSMixture) = error("NotImplemented")

"""
    components(mixture::AbstractEoSMixture)

**Required**. Iterable of corresponding eos components.
"""
components(m::AbstractEoSMixture) = error("NotImplemented")

name(m::AbstractEoSMixture) = join(map(name, components(m)), " + ")
load(::Type{T}; names, component_dbs, mix_eos_db) where {T<:AbstractEoSMixture} = error("NotImpemented")

# TODO: Perhapse, we need a (singleton) object to represent absence of buffer
"""
    thermo_buffer(mixture::AbstractEoSMixture) -> AbstractEoSThermoBuffer

Allocate buffer for intermediate calculations of mixture thermodynamic properties.

It is passed in a number of functions:
[`pressure`](@ref), [`log_c_activity!`](@ref), [`log_c_activity_wj!`](@ref)...

The method is optional, default is `nothing`.
"""
thermo_buffer(m::AbstractEoSMixture) = nothing
