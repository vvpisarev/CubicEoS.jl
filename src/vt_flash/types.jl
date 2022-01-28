abstract type AbstractVTFlashState end

"""
    value(s::AbstractVTFlashState)

Argument for optimization in VT-Flash.
"""
value(s::AbstractVTFlashState) = s.x

"""
    nmolvol(s::AbstractVTFlashState, nmolb, volumeb) -> (nmol, volume)

Moles [mol] and volume [m³] of a phase at `s`tate.
Moles `nmolb` and `volumeb` relate to base state.
"""
nmolvol(s::AbstractVTFlashState, nmolb, volumeb) = error("NotImplemented")
