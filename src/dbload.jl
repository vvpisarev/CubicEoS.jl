using CubicEoSDatabase

"""
  load(BrusilovskyEoSComponent; name, physdb)

Returns `BrusilovskyEoSComponent` by loading its parameters from database.
"""
function load(::Type{<:BrusilovskyEoSComponent}; name::AbstractString, physics_db::AbstractString, eos_db::AbstractString)
    # load physical parameters: crit temp, crit press, acentric
    comp_phys = ComponentDatabase(physics_db)
    comp_eos = ComponentDatabase(eos_db)
    
    phys = getentry(comp_phys, name)
    eos = getentry(comp_eos, name)
    return __load_brusilovsky_comp__(; phys..., eos...)
end

function __load_brusilovsky_comp__(; name::AbstractString,
                          molecular_mass,
                          number_carbons,
                          critical_temperature,
                          critical_pressure,
                          acentric_factor,
                          critical_compressibility,
                          critical_omega,
                          psi)
    __load_brusilovsky_comp__(name,
        molecular_mass,
        number_carbons,
        critical_temperature,
        critical_pressure,
        acentric_factor,
        critical_compressibility,
        critical_omega,
        psi)
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

function load(::Type{<:BrusilovskyEoSMixture};
              names,
              comp_physics_db::AbstractString,
              comp_eos_db::AbstractString,
              mix_eos_db::AbstractString)
    # load physical parameters: crit temp, crit press, acentric
    
    components = [load(BrusilovskyEoSComponent, name = name, physics_db = comp_physics_db, eos_db = comp_eos_db) for name in names]
    mix_eos = MixtureDatabase(mix_eos_db)
    corrections = getmatrix(mix_eos, names)
    return BrusilovskyEoSMixture(;
                components = components,
                corrections...)
end