"""
    molar_mass(c::AbstractEoSComponent)

Return the molar mass of a component.
"""
function molar_mass end

"""
    carbon_number(c::AbstractEoSComponent)

Return the number of carbons atoms in the hydrocarbon chain.
"""
function carbon_number end

"""
    name(c::AbstractEoSComponent)

Return the component species name.
"""
function name end

"""
    describe(c::AbstractEoSComponent)

Return `Dict` of parameters. Useful for logging.
"""
describe(x::AbstractEoSComponent) = error("NotImplemented")

"""
    load(::Type{T}; name::AbstractString, databases) where {T<:AbstractEoSComponent}

Load component `::T` by its `name` by joining the data on this species in `databases`.
"""
load(::Type{T}; name, component_dbs) where {T<:AbstractEoSComponent} = error("NotImplemented")

"""
    ncomponents(mixture::AbstractEoSMixture)

Number of components in `mixture`.
"""
function ncomponents end

@inline components(m::AbstractEoSMixture) = m.components

@inline ncomponents(m::AbstractEoSMixture) = length(components(m))
@inline Base.length(m::AbstractEoSMixture) = length(components(m))

# Base.length(m::AbstractEoSMixture) = m.number_of_components
Base.show(io::IO, x::AbstractEoSMixture) = print(io, "$(typeof(x))($(name(x)))")

name(m::AbstractEoSMixture) = join(map(name, components(m)), " + ")
describe(m::AbstractEoSMixture) = Dict{String,Any}("noparameters" => NaN)

load(::Type{T}; names, component_dbs, mix_eos_db) where {T<:AbstractEoSMixture} = error("NotImpemented")
