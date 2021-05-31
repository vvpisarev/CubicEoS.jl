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
    0.01 0
]

c1c4mix = BrusilovskyEoSMixture(
    components = [ch4, c4h10],
    constant = c1c4cij,
    linear = zeros(2,2),
    quadratic = zeros(2,2)
)

nmol = [500., 500]

_log_a = zeros(2)
_ai = zeros(2)
_aij = zeros(2, 2)

p = pressure(c1c4mix, nmol, 2.0, 8.314 * 300; ai = _ai, aij = _aij)
println("Pressure: ", p, " Pa")

function log_a(mix, nmol, vol, RT)
    return CubicEoS.log_activity(mix, nmol, vol, RT; ai = _ai, aij = _aij)
end

chempot(nmol) = 8.314 .* 300 .* (log.(nmol) .+ log_a(c1c4mix, nmol, 2.0, 8.314*300))

function helmholtz(nmol, vol)
    f = sum(nmol .* chempot(nmol))
    f -= vol * pressure(c1c4mix, nmol, vol, 8.314 * 300; ai = _ai, aij = _aij)
    return f
end

println("Helmholtz energy: ", helmholtz(nmol, 2.0))

nmol1 = nmol .+ (1, 0)
println("Approx. chemical potential: ", helmholtz(nmol1, 2.0) - helmholtz(nmol, 2.0))

println("Chemical potential: ", chempot(nmol))

(log_a(c1c4mix, nmol1, 2.0, 8.314*300) - log_a(c1c4mix, nmol, 2.0, 8.314*300)) |> println

CubicEoS.log_activity_wj(
    c1c4mix,
    nmol,
    2.0,
    8.314*300
)
