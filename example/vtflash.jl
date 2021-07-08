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

V = 1e-6
χ = [0.547413, 0.452587]

for T in 300:10:450
    for Ση in 100:100:15000
        RT = T * CubicEoS.GAS_CONSTANT_SI
        N = (Ση * V) .* χ
        converged = false
        singlephase = false
        try
            state = vt_flash(C₁C₅, N, V, RT)
            converged = state.converged
            singlephase = state.singlephase
        catch e
            @warn "VTFlash not converged" T Ση e
        end
        print(join([T, Ση], '\t'))
        print('\t', join([converged, singlephase], '\t'))
        println()
    end
end
