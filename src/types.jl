"""
    AbstractEoSComponent

Abstract type representing component.
"""
abstract type AbstractEoSComponent end

"""
    AbstractEoSMixture{T}

Abstract type representing mixture.

Parameter `T` should represents `<:Real` subtype to store fractional values like
molar mass, critical temperature and pressure, etc.
However, no type-constraints added by `CubicEoS` to support generics (e.g. dual numbers).
"""
abstract type AbstractEoSMixture{T} end

"""
    AbstractEoSThermoBuffer

Abstract type representing buffer for intermediate calculations.

See [`thermo_buffer`](@ref).
"""
abstract type AbstractEoSThermoBuffer end


"""
    OptimStats(converged, iters, calls)

Results of optimization.

# Getters

- `converged(stats)::Bool`: if `true`, optimization succesfully finished and converged;
- `iters(stats)::Int`: number of optimization steps (iterations) done;
- `calls(stats)::Int`: number of function calls done.
"""
struct OptimStats
    converged::Bool
    iters::Int
    calls::Int
end

converged(x::OptimStats) = x.converged
iters(x::OptimStats) = x.iters
calls(x::OptimStats) = x.calls
