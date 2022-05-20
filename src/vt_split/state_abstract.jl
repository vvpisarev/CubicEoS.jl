"""
    value(s::AbstractVTSplitState)

Argument for optimization in VT-Flash.
"""
value(s::AbstractVTSplitState) = s.x

"Constructor of a VTFlashState from concentration and saturation."
(::Type{<:AbstractVTSplitState})(concentration, saturation, nmolbase, volumebase) = error("NotImplemented")

"""
    nmolvol!(nmol, state::AbstractVTSplitState, nmolbase, volumebase) -> (nmol, volume)

Destructive option of [`nmolvol`](@ref).
"""
nmolvol!(nmol, s::AbstractVTSplitState, nmolbase, volumebase) = error("NotImplemented")

"""
    nmolvol(state::AbstractVTSplitState, nmolbase, volumebase) -> (nmol, volume)

Converts a `state` variables into moles `nmol` and `volume`
from `state` of a corresponding phase.

See also [`nmolvol!`](@ref).
"""
function nmolvol(s::AbstractVTSplitState, nmolb::V, volumeb::T) where {V, T}
    nmol = similar(nmolb, Float64)  # FIXME: type inference
    return nmolvol!(nmol, s, nmolb, volumeb)
end

"Vector of upper limits for physical constraint for given state variables."
physical_constrain_step_uplims(::Type{<:AbstractVTSplitState}, nmolbase, volumebase) = error("NotImplemented")

# Creates closure of physical constrains `0 < x + α d < xuplims`.
function physical_constrain_step_closure(
    ::Type{<:AbstractVTSplitState},
    xuplims::AbstractVector,
)
    function clsr(
        ::Type{<:AbstractVTSplitState},
        x::AbstractVector,
        direction::AbstractVector,
    )
        αmax = Inf
        for (i, (xi, di, ui)) in enumerate(zip(x, direction, xuplims))
            if iszero(di) && !(0 < xi < ui)
                throw(ConstrainStepZeroDirectionError(i, xi))
            end
            αmax = di < 0 ? min(αmax, -xi/di) : min(αmax, (ui - xi)/di)
        end

        αlow = -Inf
        for (xi, di, ui) in zip(x, direction, xuplims)
            # Zero direction is checked above
            αlow = di > 0 ? max(αlow, -xi/di) : max(αlow, (ui - xi)/di)
        end
        return αlow, αmax
    end
    return clsr
end
