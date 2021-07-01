#=
VT-flash algorithm
=#

struct VTFlashResult{T}
    converged::Bool
    singlephase::Bool
    RT::T
    nmol_1::Vector{T}
    V_1::T
    nmol_2::Vector{T}
    V_2::T
end

VTFlashResult{T}(; converged, singlephase, RT, nmol_1, V_1, nmol_2, V_2) where {T} =
VTFlashResult{T}(converged, singlephase, RT, nmol_1, V_1, nmol_2, V_2)

function vt_flash(
    mix::BrusilovskyEoSMixture{T},
    nmol::AbstractVector,
    volume::Real,
    RT::Real,
) where {T}
    # run vt-stability to find out whether a state single phase or not
    singlephase, η_base = vt_stability(mix, nmol, volume, RT)
    if singlephase
        return VTFlashResult{T}(
            converged=true,
            singlephase=true,
            RT=RT,
            nmol_1=nmol,
            V_1=volume,
            nmol_2=similar(nmol),
            V_2=0
        )
    end

    # two-phase state case
    # create closures for helmoltz energy, its gradient and constrain step

    # find initial vector for optimizer

    # run optimizer

    # sort phases into foo_1 for gas, foo_2 for liquid

    # return result

    return nothing
end
