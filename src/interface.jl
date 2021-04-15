export molar_mass, carbon_number, name, describe, load, components
"""
    molar_mass(c::AbstractEoSComponent)

Return the molar mass of a component.
"""
function molar_mass end

"""
    molar_mass(c::AbstractEoSComponent)

Return the molar mass of a component.
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
function describe end

"""
    load(::Type{T}; name::AbstractString, physdb::AbstractString) where {T<:AbstractEoSComponent}

Load component `::T` by its `name` and `physdb` from database.
"""
function load end


Base.length(m::AbstractEoSMixture) = m.number_of_components
Base.show(io::IO, x::AbstractEoSMixture) = print(io, "$(typeof(x))($(name(x)))")

components(m::AbstractEoSMixture) = m.components
name(m::AbstractEoSMixture) = join(map(name, components(m)), " + ")
describe(m::AbstractEoSMixture) = Dict{String,Any}("noparameters" => NaN)
load(::Type{T}; names, physdb) where {T<:AbstractEoSMixture} = error("NotImpemented")

function describe(x::BrusilovskyEoSComponent)
    return Dict{String,Any}(
        "data structure" => repr(x),
        "name" => name(x),
        "critical pressure [Pa]" => x.Pc,
        "critical temperature [K]" => x.Tc,
        "pitzer acentric factor" => x.acentric_factor,
        "molar mass [kg mol⁻¹]" => x.molar_mass,
        "number of carbons atoms" => x.carbon_number,
        "eos" => "brusilovsky",
        "eos param: ac [?]" => x.ac,
        "eos param: b [?]"  => x.b,
        "eos param: c [?]"  => x.c,
        "eos param: d [?]"  => x.d,
    )
end