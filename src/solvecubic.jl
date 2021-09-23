"""
    solve_cubic(a::Real, b::Real, c::Real, d::Real)

Finds real roots of a cubic equation
```math
a x^3 + b x^2 + c x + d = 0
```

Return Tuple{Float64, Float64, Float64} where first roots are real,
complex roots are defined to `NaN`

Reference: J. F. Blinn, "How to Solve a Cubic Equation, Part 5: Back to Numerics," in IEEE Computer Graphics and Applications, vol. 27, no. 3, pp. 78-89, May-June 2007.
doi: 10.1109/MCG.2007.60
"""
function solve_cubic(a::Real, b::Real, c::Real, d::Real)
    # convert to Ax³ + 3Bx² + 3Cx + D = 0
    A, B, C, D = promote(float(a), float(b) / 3, float(c) / 3, float(d))
    δ₁ = A * C - B * B
    δ₂ = A * D - B * C
    δ₃ = B * D - C * C
    d13 = δ₁ * δ₃
    d22 = δ₂ * δ₂
    Δ = 4 * d13 - d22
    if abs(Δ) < eps(d22)
        Δ = zero(Δ)
    end
    nanvalue = zero(A) / zero(A)

    if Δ < 0
        At, Cb, Db = zero.((A, C, D)) # A-tilde, C-bar, D-bar
        if B^3 * D >= A * C^3
            At, Cb, Db = A, δ₁, -2 * B * δ₁ + A * δ₂
        else
            At, Cb, Db = D, δ₃, -D * δ₂ + 2 * C * δ₃
        end
        T₀ = -copysign(At, Db) * sqrt(-Δ)
        T₁ = -Db + T₀
        p = cbrt(T₁ / 2)
        q = T₁ == T₀ ? -p : -Cb / p
        x₁ = Cb <= 0 ? p + q : -Db / (p^2 + q^2 + Cb)
        x, w = B^3 * D >= A * C^3 ? (x₁ - B, A) : (-D, x₁ + C)
        return (x / w, nanvalue, nanvalue)
    else
        δ₁ == δ₂ == δ₃ == 0 && return (-B/A, -B/A, -B/A)
        sΔ = sqrt(Δ)
        θA, θD = (atan(A*sΔ, 2*B*δ₁ - A*δ₂), atan(D*sΔ, D*δ₂ - 2*C*δ₃)) ./ 3 .|> abs
        sCA, sCD = sqrt(-δ₁), sqrt(-δ₃)
        x₁A, x₁D = 2*sCA*cos(θA), 2*sCD*cos(θD)
        x₃A, x₃D = -sCA*(cos(θA)+sqrt(3)*sin(θA)), -sCD*(cos(θD)+sqrt(3)*sin(θD))
        xlt = (x₁A+x₃A > 2*B) ? x₁A : x₃A
        xst = (x₁D+x₃D < 2*C) ? x₁D : x₃D
        xl, wl = xlt - B, A
        xs, ws = -D, xst + C
        Δ == 0 && return (xs / ws, xl / wl, nanvalue)
        E = wl * ws
        F = -xl * ws - wl * xs
        G = xl * xs
        xm, wm = C * F - B * G, C * E - B * F
        return (xs / ws, xm / wm, xl / wl)
    end
end

"""
    solve_cubic!(roots!::AbstarctVector, a, b, c, d) -> roots!

See [`solve_cubic`](@ref). Fill first three items of `roots!`.
"""
function solve_cubic!(roots!::AbstractVector, a::Real, b::Real, c::Real, d::Real)
    roots![1:3] .= solve_cubic(a, b, c, d)
    return roots!
end
