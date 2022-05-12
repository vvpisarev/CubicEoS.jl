CubicEoS.molar_mass(x::BrusilovskyEoSComponent) = x.molar_mass
CubicEoS.name(x::BrusilovskyEoSComponent) = x.name
CubicEoS.carbon_number(x::BrusilovskyEoSComponent) = x.carbon_number

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
