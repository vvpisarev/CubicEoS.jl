using CubicEoS
using CubicEoSDatabase: Data

C₁C₅ = load(
    BrusilovskyEoSMixture,
    names = ("methane", "n-pentane"),
    component_dbs = (Data.martinez(), Data.brusilovsky_comp()),
    mix_eos_db = Data.brusilovsky_mix()
)

# [0.547413, 0.452587] mole fractions from (Mikyska, 2013, Example 2)
# V = 1e-6 from old vtflash
# (5e3 * V) is (5 kmol m⁻³ * V) just an ordinary value
# T = 370 Kelvins at 5 kmol m⁻³ should correspond to two-phase
# T = 440 Kelvins at 5 kmol m⁻³ should correspond to single-phase

V = 1#1e-6
N = (5000 * V) .* [0.547413, 0.452587]
RT = 300 * CubicEoS.GAS_CONSTANT_SI

issinglephase, tries = vt_stability(C₁C₅, N, V, RT)
dump(tries)

dump(CubicEoS.__vt_flash_init_conc_choose(tries))