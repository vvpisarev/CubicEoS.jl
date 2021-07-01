using CubicEoS
using CubicEoSDatabase: Data


c1c5mix = load(
    BrusilovskyEoSMixture,
    names = ("methane", "n-pentane"),
    component_dbs = (Data.martinez(), Data.brusilovsky_comp()),
    mix_eos_db = Data.brusilovsky_mix()
)

"""
Prints (Ση, T)-diagram (variables used in (Mikyska, 2012)) of single phase stability region.

Uses VT-stability solver.

The mixture of composition `χ_base` is loaded by `names` of its components.
Then calculation performed for each (T ∈ `T_range`, Ση_base ∈ `Ση_range`) point.

Warns through stderr if an exception happens.
"""
function vt_stability_η_temp(io::IO, χ_base, Ση_range, T_range)
    mixture = c1c5mix

    let
        prefix = "# "
        println(io, prefix, "VT-stability calculation of (Ση, T)-diagram")
        println(io, prefix, "Mixture: $(name(mixture))")
        println(io, prefix, "Composition (in molar parts) of base phase χ_base: ")
        for (c, χ) in zip(components(mixture), χ_base)
            println(io, prefix, "  $(name(c))\t\t$χ")
        end
        println(io, prefix, "Range of overall concentration Ση (mol m⁻³):")
        println(io, prefix, "  from: ", first(Ση_range))
        println(io, prefix, "    to: ", last(Ση_range))
        println(io, prefix, "Range of temperatures T (K):")
        println(io, prefix, "  from: ", first(T_range))
        println(io, prefix, "    to: ", last(T_range))
        
        columns = (
            ("T", "temperature (К)"),
            ("Ση", "overall concentration of mixture (mol m⁻³)"),
            # generator of tuples ("η_base_i", "description of η_base_i")
            (("η_base_$i", "concentration (mol m⁻³) of $(name(c)) in base phase") for (i, c) in enumerate(components(mixture)))...,
            ("is_fail", "if true, an error raised during calculation"),
            ("is_stable", "if true, base phase is stable (single phase state)"),
        )
        println(io, prefix, "Data representation")
        println(io, prefix, "  <column index>. <symbol>: <description>")
        for (i, col) in enumerate(columns)
            println(io, prefix, "  $i. ", col[1], ": ", col[2])
        end
    end

    sep = '\t'

    η_base = fill(NaN, length(mixture))

    for T in T_range
        for Ση in Ση_range
            η_base .= Ση .* χ_base
            RT = T * CubicEoS.GAS_CONSTANT_SI

            is_fail = false
            is_stable = false
            
            try
                is_stable, = vt_stability(mixture, η_base, 1.0, RT)
            catch e
                @warn "Exception occurs" exception=e η_base=repr(η_base) RT=RT
                is_fail = true
                is_stable = false
            finally
                println(io, join([T, Ση, η_base..., repr(is_fail), repr(is_stable)], sep))
            end
        end
    end
end

vt_stability_η_temp(χ_base, Ση_range, T_range) = vt_stability_η_temp(stdout, χ_base, Ση_range, T_range)

function vt_stability_η_temp(fname::AbstractString, χ_base, Ση_range, T_range)
    open(fname, "w") do io
        vt_stability_η_temp(io, χ_base, Ση_range, T_range)
    end
end

vt_stability_η_temp([0.547413, 0.452587], 100:10:15000, 300:1:450)
