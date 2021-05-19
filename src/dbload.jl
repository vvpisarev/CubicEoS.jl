using CubicEoSDatabase

"""
    load(BrusilovskyEoSComponent; name::AbstractString[, physics_db::ComponentDatabase, eos_db::ComponentDatabase])

Returns `BrusilovskyEoSComponent` named `name` by loading its parameters from

- database of physical properties `physics_db` (default is `CubicEoSDatabase.Data.martinez()`);
- database of eos properties `eos_db` (default is `CubicEoSDatabase.Data.brusilovsky_comp()`).
"""
function load(::Type{<:BrusilovskyEoSComponent};
        name::AbstractString,
        physics_db::ComponentDatabase=Data.martinez(),
        eos_db::ComponentDatabase=Data.brusilovsky_comp()
    )
    phys = getentry(physics_db, name)
    eos = getentry(eos_db, name)
    return __load_brusilovsky_comp__(; phys..., eos...)
end

function __load_brusilovsky_comp__(; name::AbstractString,
                          molecular_mass,
                          number_carbons,
                          critical_temperature,
                          critical_pressure,
                          acentric_factor,
                          eos_critical_compressibility,
                          eos_critical_omega,
                          eos_psi,
                          extra_kw...  # unneeded keywords to construct component
    )
    __load_brusilovsky_comp__(name,
        molecular_mass,
        number_carbons,
        critical_temperature,
        critical_pressure,
        acentric_factor,
        eos_critical_compressibility,
        eos_critical_omega,
        eos_psi)
end

function __load_brusilovsky_comp__(name::AbstractString,
            molecular_mass,
            number_carbons,
            critical_temperature,
            critical_pressure,
            acentric,
            critical_compressibility::AbstractString,
            critical_omega::AbstractString,
            psi::AbstractString)

    Omegac  = 0.75001                       # Table (4.14)
    Zc      = 0.3357 - 0.0294 * acentric    # 4.203
    if acentric < 0.4489
        Psi = 1.050 + ( 0.105 + 0.482 * acentric ) * acentric   # 4.204
    else
        Psi = 0.429 + ( 1.004 + 1.561 * acentric ) * acentric   # 4.204
    end
    return BrusilovskyEoSComponent(
        name = name,
        critical_pressure = critical_pressure,
        critical_temperature = critical_temperature,
        acentric_factor = acentric,
        Omegac = Omegac,
        Zc = Zc,
        Psi = Psi,
        molar_mass = molecular_mass,
        carbon_number = number_carbons
    )
end

function __load_brusilovsky_comp__(name::AbstractString,
            molecular_mass,
            number_carbons,
            critical_temperature,
            critical_pressure,
            acentric,
            critical_compressibility,
            critical_omega,
            psi)

    return BrusilovskyEoSComponent(
                name = name,
                critical_pressure = critical_pressure,
                critical_temperature = critical_temperature,
                acentric_factor = acentric,
                Omegac = critical_omega,
                Zc = critical_compressibility,
                Psi = psi,
                molar_mass = molecular_mass,
                carbon_number = number_carbons
            )
end

"""
    load(BrusilovskyEoSMixture; names[, comp_physics_db::ComponentDatabase, comp_eos_db::ComponentDatabase, mix_eos_db::MixtureDatabase])

Returns `BrusilovskyEoSMixture` of `names` by loading parameters from

- database of physical properties `physics_db`. Default is `CubicEoSDatabase.Data.martinez()`;
- database of eos properties `eos_db`. Default is `CubicEoSDatabase.Data.brusilovsky_comp()`;
- database of eos binary interaction parameters `mix_eos_db`. Default is `CubicEoSDatabase.Data.brusilovsky_mix()`.
"""
function load(::Type{<:BrusilovskyEoSMixture};
              names,
              comp_physics_db::ComponentDatabase=Data.martinez(),
              comp_eos_db::ComponentDatabase=Data.brusilovsky_comp(),
              mix_eos_db::MixtureDatabase=Data.brusilovsky_mix())
    
    components = [load(BrusilovskyEoSComponent, name = name, physics_db = comp_physics_db, eos_db = comp_eos_db) for name in names]
    corrections = getmatrix(mix_eos_db, names)
    return BrusilovskyEoSMixture(;
                components = components,
                corrections...)
end