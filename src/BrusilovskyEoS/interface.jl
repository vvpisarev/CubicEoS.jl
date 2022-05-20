CubicEoS.molar_mass(x::BrusilovskyEoSComponent) = x.molar_mass
CubicEoS.name(x::BrusilovskyEoSComponent) = x.name
CubicEoS.carbon_number(x::BrusilovskyEoSComponent) = x.carbon_number

CubicEoS.components(x::BrusilovskyEoSMixture) = x.components
CubicEoS.ncomponents(x::BrusilovskyEoSMixture) = length(components(x))
Base.show(io::IO, x::BrusilovskyEoSMixture) = print(io, "BrusilovskyEoSMixture($(CubicEoS.name(x)))")
CubicEoS.thermo_buffer(mix::BrusilovskyEoSMixture) = BrusilovskyThermoBuffer(mix)
