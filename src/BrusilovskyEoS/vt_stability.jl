function eos_constrain_step(
    StateVariables::Type{CubicEoS.VTStabilityIdealIdentityState},
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
