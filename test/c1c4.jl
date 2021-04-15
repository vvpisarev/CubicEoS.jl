ch4 = BrusilovskyEoSComponent(
    name = "CH4",
    critical_pressure = 4.5992e6,
    critical_temperature = 190.564,
    acentric_factor = 0.01142,
    Omegac = 0.7563,
    Zc = 0.33294,
    Psi = 0.37447,
    carbon_number = 1,
    molar_mass = 12.011 + 4 * 1.008
)

c4h10 = BrusilovskyEoSComponent(
    name = "C4H10",
    critical_pressure = 3.7960e6,
    critical_temperature = 425.125,
    acentric_factor = 0.201,
    Omegac = 0.76921,
    Zc = 0.31232,
    Psi = 0.57594,
    carbon_number = 4,
    molar_mass = 12.011 * 4 + 10 * 1.008
)

c1c4cij = [
    0 0.01
    0.01 0]

c1c4mix = BrusilovskyEoSMixture(components = [ch4, c4h10],
                                constant = c1c4cij,
                                linear = zeros(2,2),
                                quadratic = zeros(2,2))

c1c4_chem_mix = CubicEoS.ChemPotentialMixture(c1c4mix)

nmol = [500., 500]

p = pressure(c1c4mix, nmol, 2.0, 8.314 * 300)
println("Pressure: ", p, " Pa")

log_a = zeros(2)
ai = zeros(2)
aij = zeros(2, 2)

chempot(nmol) = 8.314 .* 300 .* (log.(nmol) .+ 
                                 CubicEoS.log_activity!(c1c4_chem_mix, 
                                                        nmol, 
                                                        2.0, 
                                                        8.314*300))
helmholtz(nmol) = sum(nmol .* chempot(nmol)) - CubicEoS.pressure!(c1c4_chem_mix, nmol, 2.0, 8.314 * 300)
println("Helmholtz energy: ", helmholtz(nmol))

nmol1 = nmol .+ (1, 0)
println("Approx. chemical potential: ", helmholtz(nmol1) - helmholtz(nmol))

println("Chemical potential: ", chempot(nmol))

(CubicEoS.log_activity(
c1c4mix, 
nmol1, 
2.0, 
8.314*300) - CubicEoS.log_activity(
    c1c4mix, 
    nmol, 
    2.0, 
    8.314*300)) |> println

CubicEoS.log_activity_wj!(
c1c4_chem_mix, 
nmol, 
2.0, 
8.314*300)