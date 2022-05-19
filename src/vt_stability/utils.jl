"EoS constraint with no bounds on magnitude of step for vt stability."
unconstrained_eos_step(::Type{<:AbstractVTStabilityState}, x, d) = (-Inf, Inf)
