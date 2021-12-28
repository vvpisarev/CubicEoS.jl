using CubicEoS
# using CubicEoSDatabase: Data
using LinearAlgebra


function vt_flash_helmoltz_diff_saturation_tries(; mix, nmol, volume, RT)
    issinglephase, vt_stab_tries = vt_stability(mix, nmol, volume, RT)
    @assert !issinglephase

    function get_state(conc, sat)
        return [(conc * (volume * sat) ./ nmol)..., sat]
    end

    grad = fill(NaN, ncomponents(mix) + 1)
    for sat₁ in 0.01:0.01:0.99
        ΔA_vec = fill(NaN, length(vt_stab_tries))
        error_vec = fill(false, length(vt_stab_tries))
        for (j, vt_try) in enumerate(vt_stab_tries)
            state = get_state(vt_try.concentration, sat₁)
            ΔA, error = NaN, false
            try
                ΔA, grad = helmholtz_diff!(state, grad)
            catch e
                @warn "helmholtz_diff! error" sat₁ e
                ΔA, error = NaN, true
            end
            ΔA_vec[j] = ΔA
            error_vec[j] = error
        end
        print(sat₁)
        print('\t', join(error_vec, '\t'))
        print('\t', join(ΔA_vec, '\t'))
        println()
    end
    return nothing
end

function vt_flash_helmoltz_diff_saturation(
    ;
    mix,
    nmol,
    volume,
    RT,
    conc₁,
    helmholtz_diff!::Function,
    satinit,
)
    state = fill(NaN, ncomponents(mix) + 1)
    grad = similar(state)

    function upd_state(conc₁, sat₁)
        N₁ = conc₁ * (volume * sat₁)
        state[1:end-1] .= N₁ ./ nmol
        state[end] = sat₁
        return nothing
    end

    sat₁ = satinit
    scale = 1
    scale_step = 0.9
    upd_state(conc₁, sat₁)
    for step in 1:50
        dir = 1 * [(conc₁ * (volume * sat₁) ./ nmol)..., 1]
        dir ./= norm(dir, 2)
        scalemax = constrain_step(state, dir) / 0.9
        satmax = state[end] + scalemax * dir[end]

        sat₁ = satinit * scale
        upd_state(conc₁, sat₁)
        ΔA, error = NaN, false
        try
            ΔA, grad = helmholtz_diff!(state, grad)
        catch e
            @warn "helmholtz_diff! error" sat₁ e
            ΔA, error = NaN, true
        end
        scale *= scale_step

        print(join([sat₁, satmax], '\t'))
        print('\t', error)
        print('\t', join([ΔA, norm(grad)], '\t'))
        println()
    end
    # upd_state(conc₁, 0.8)
    # ΔA, grad = helmholtz_diff!(state, grad)
    return nothing
end

mix = load(
    BrusilovskyEoSMixture;
    names=("methane", "n-pentane"),
    # component_dbs=(Data.martinez(), Data.brusilovsky_comp()),
    # mix_eos_db=Data.brusilovsky_mix()
)
# dump(mix)
# # exit()
# [0.547413, 0.452587] mole fractions from (Mikyska, 2013, Example 2)
# V = 1e-6 from old vtflash
# (5e3 * V) is (5 kmol m⁻³ * V) just an ordinary value
# T = 370 Kelvins at 5 kmol m⁻³ should correspond to two-phase
# T = 440 Kelvins at 5 kmol m⁻³ should correspond to single-phase

V = 1e-6
N = (5000 * V) .* [0.547413, 0.452587]
RT = 370 * CubicEoS.GAS_CONSTANT_SI

# closures
constrain_step, helmholtz_diff_grad!, helmholtz_diff! = CubicEoS.vt_flash_closures(mix, N, V, RT)

# vtstability
issinglephase, vt_stab_tries = vt_stability(mix, N, V, RT)
conc₁ = CubicEoS.__vt_flash_init_conc_choose(vt_stab_tries)
println(stderr, "Is single phase: ", issinglephase)
println(stderr, "Concentration from vt_stability ", conc₁)
# for vttry in vt_stab_tries
#     dump(vttry)
# end

# initial state debug
# state = fill(NaN, ncomponents(mix) + 1)
# found_initial = CubicEoS.vt_flash_initial_state!(
#     state, N, V, conc₁, helmholtz_diff!, constrain_step;
#     sat₁max=0.1,
#     steps=100,
#     step_scale=0.8,
#     helmholtz_thresh=-1e-5,
# )
# println(stderr, "vt_flash_initial_state", " state is found: ", found_initial)
# println(helmholtz_diff!(state, similar(state)))

# ΔA from saturation and initial concentration
# vt_flash_helmoltz_diff_saturation(
#     ;
#     mix=mix,
#     nmol=N,
#     volume=V,
#     RT=RT,
#     conc₁=conc₁,
#     helmholtz_diff! =helmholtz_diff!,
#     satinit=1,
# )

# vt_flash_helmoltz_diff_saturation_tries(mix=mix, nmol=N, RT=RT, volume=V)

# ∇P = fill(NaN, ncomponents(mix) + 1)
# CubicEoS.__vt_flash_pressure_gradient!(∇P, mix, N, V, RT)
# dump(∇P)

# ℍ = fill(NaN, (ncomponents(mix) + 1, ncomponents(mix) + 1))
# # state = fill(NaN, size(ℍ)[1])
# CubicEoS.__vt_flash_hessian!(ℍ, state, mix, N, V, RT)
# dump(ℍ)

state = CubicEoS.vt_flash_newton(mix, N, V, RT)
dump(state)
