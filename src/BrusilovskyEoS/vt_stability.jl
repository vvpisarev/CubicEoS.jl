function eos_vt_stability_constrain_step(
    ::Type{<:CubicEoS.VTStabilityIdealIdentityState},
    mix::BrusilovskyEoSMixture,
)
    covolumes = components(mix).b

    function clsr(StateVariables, trialx, direction)
        # Size of trial phase
        A = sum(d*d*b for (d, b) in zip(direction, covolumes))
        B = 2 * sum(x*d*b for (x, d, b) in zip(trialx, direction, covolumes))
        C = -4 + sum(x*x*b for (x, b) in zip(trialx, covolumes))

        A < 0 && throw(CubicEoS.ConstrainStepError("parabolic condition wrong"))
        αcov = CubicEoS.solve_quadratic(A, B, C)

        # If there are one or zero roots, covolume inequality constraint has no solutions.
        any(isnan, αcov) && throw(CubicEoS.ConstrainStepError("covolume constraint can't be resolved"))

        αlow, αmax = extrema(αcov)

        return αlow, αmax
    end

    return clsr
end

# Condition: dot(x + α d, b) < 1,
# where x is concentration, d is direction and b is covolumes.
function eos_vt_stability_constrain_step(
    ::Type{<:CubicEoS.VTStabilityPhysicalState},
    mix::BrusilovskyEoSMixture,
)
    covolumes = components(mix).b

    function clsr(StateVariables, x, direction)
        dirdotcov = dot(direction, covolumes)
        if iszero(dirdotcov)
            (1 - dot(x, covolumes)) < 0 && throw(CubicEoS.ConstrainStepError("covolume constraint can't be resolved"))
            # Correct dot(x, covolumes)
            return -Inf, Inf
        end

        αmax = dirdotcov > 0 ? (1 - dot(x, covolumes)) / dirdotcov : Inf
        αlow = dirdotcov < 0 ? (1 - dot(x, covolumes)) / dirdotcov : -Inf
        return αlow, αmax
    end
end
