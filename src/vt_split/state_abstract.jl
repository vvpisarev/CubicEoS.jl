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

# TODO: refactor names of `physical_constrain_step_*lims`: low/high, down/up, bottom/top?..

"""
    CubicEoS.physical_constrain_step_lowlims(::Type{<:AbstractVTSplitState}, nmolbase, volumebase)
    -> Vector of size size(nmolbase, 1) + 1

Vector of lower limits for physical constraint for given state variables.

By default, consists of zeros. This is suitable, for example, when state is represented by moles and volume.

See also [`CubicEoS.physical_constrain_step_uplims`](@ref).
"""
physical_constrain_step_lowlims(::Type{<:AbstractVTSplitState}, nmolbase, volumebase) = zeros(size(nmolbase, 1) + 1)

"""
    CubicEoS.physical_constrain_step_uplims(::Type{<:AbstractVTSplitState}, nmolbase, volumebase)

**Must be implemented**. Vector of upper limits for physical constraint for given state variables.

For example, when state is represented by moles and volume, the uplims are `[nmolbase; volumebase]`.

See also [`CubicEoS.physical_constrain_step_lowlims`](@ref).
"""
physical_constrain_step_uplims(::Type{<:AbstractVTSplitState}, nmolbase, volumebase) = error("NotImplemented: CubicEoS.physical_constrain_step_uplims")


# Creates closure of physical constraints `xlowlims < x + α d < xuplims`.
# The closure defines bounds for step magnitude.
function physical_constrain_step_closure(
    ::Type{<:AbstractVTSplitState},
    xlowlims::AbstractVector,
    xuplims::AbstractVector,
)
    function clsr(
        ::Type{<:AbstractVTSplitState},
        x::AbstractVector,
        direction::AbstractVector,
    )
        # Solving of low[i] < x[i] + α * d[i] < high[i] for α and all i-s.

        αmax, αlow = Inf, -Inf
        for (i, (xi, di, li, ui)) in enumerate(zip(x, direction, xlowlims, xuplims))
            if iszero(di) && !(li < xi < ui)
                # Zero direction, but current value is outside bounds.
                throw(ConstrainStepZeroDirectionError(i, xi))
            end

            lbound, ubound = (li - xi)/di, (ui - xi)/di

            αlow_candidate = di > 0 ? lbound : ubound
            αlow = max(αlow, αlow_candidate)

            αmax_candidate = di < 0 ? lbound : ubound
            αmax = min(αmax, αmax_candidate)
        end
        return αlow, αmax
    end
    return clsr
end
