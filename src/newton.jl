using LinearAlgebra

"""
Result of Newton's iterations.
Contain `converged::Bool` flag, the minimum `argument::Vector`,
number of `iterations::Integer` and `calls::Integer`.
"""
struct NewtonResult{T}
    converged::Bool
    argument::Vector{T}
    iterations::Int
    calls::Int

    NewtonResult(conv, x, iters=-1, calls=-1) = new{Float64}(conv, copy(x), iters, calls)
end

"""
    newton(f, grad!, hess!, x₀[; gtol=1e-6, convcond, maxiter=200, constrain_step])

Find (un)constrained minimum of `f(x)::Real` using Newton's solver with line search
and hessian modified by Cholesky factorization.
The optimization starts at `x₀` performing no more than `maxiter` Newton's steps.

New point `x ← xpre + αmax * d` is constrained in Newton's `d`irection by
`constrain_step(x, d) -> α::Real`.
By default, the constrains are disabled.

Return [`NewtonResult`](@ref) object.

# Arguments

The target function and initialization of optimization is defined by following arguments

- `f(x) -> Real`: function to optimize;
- `grad!(g, x) -> g`: must update and return `g`radient of `f` at `x`;
- `hess!(H, x) -> H`: must update and return `H`essian matrix of `f` at `x`.
- `x₀::AbstractVector`: initial argument for optimization.

# Optional arguments

## Convergence criteria

By default, the optimization considered converged when `norm(∇(x), 2) ≤ gtol`.
Generally, the convergence criterion is defined by `convcond`.

- `gtol::Real=1e-6`: (default criterion) 2-norm of gradient when iterations
    are considered converged;
- `convcond::Function` with signature `(x, xpre, y, ypre, g)->Bool`:
    if returns `true` for given arguments `x` and `xpre`, function values `y=f(x)`,
    `ypre=f(xpre)` and gradient `g=∇(x)`, the optimization is stopped.

Example (default criterion): `stopcond=(x, xpre, y, ypre, g)->norm(g, 2) ≤ gtol`.

## Defining constrains

- `constrain_step::Function=(x, d)->Inf`:
    given an argument `x` and `d`irection, returns maximum allowed scalar `α` in form of
    `xnew = x + α*d`. Default is unconstrained optimization (`α = Inf`).

## Forcing stop of optimization

- `maxiter::Integer=200`: maximum allowed Newton's steps.
"""
function newton(
    f::Function,
    ∇f!::Function,
    H!::Function,
    x_::AbstractVector;
    gtol::Real=1e-6,
    maxiter::Integer=200,
    constrain_step::Function=(x, d)->Inf,
    convcond::Function=(x, xpre, y, ypre, g)->norm(g, 2) ≤ gtol,
)
    x = float.(x_)  # copy x for internal use

    ∇ = similar(x)
    δx = similar(x)
    hess_full = Matrix{eltype(x)}(undef, size(∇, 1), size(∇, 1))
    vec = similar(x)
    xpre = similar(x)
    totfcalls = 0


    """
    Returns closure for Downhill.strong_backtracking!
    Vector `v` for determine types.
    """
    function fdfclosure(v::AbstractVector)
        vecx = similar(v)
        vecg = similar(v)

        function fdf(x::AbstractVector, α::Real, d::AbstractVector)
            vecx .= x .+ α .* d
            return f(vecx), ∇f!(vecg, vecx)
        end

        return fdf
    end

    fdf = fdfclosure(x)

    # Initializing, main cycle is postconditioned
    y = f(x)
    totfcalls += 1
    ∇ = ∇f!(∇, x)

    for i in 1:maxiter
        # descent direction `δx`

        hess_full = H!(hess_full, x)
        hess = Downhill.mcholesky!(hess_full)

        vec .= ∇
        ldiv!(hess, vec)  # `hess \ ∇`, result in `vec`
        @. δx = -vec

        # choosing step in `δx`
        αmax = constrain_step(x, δx)
        α = Downhill.strong_backtracking!(fdf, x, δx;
            α=1.0,  # the function takes care of initial `α` and `αmax`
            αmax=αmax,
            σ=0.9,
        )
        fcalls = 0  # by current implementation they are unknown :c
        totfcalls += fcalls

        if α ≤ 0
            @error "newton: linesearch failed" α
            return NewtonResult(false, x, i, totfcalls)
        end

        # update argument `x`, function value and gradient in `x`
        xpre .= x
        ypre = y

        @. x += α * δx
        y = f(x)
        ∇ = ∇f!(∇, x)

        @debug "newton" i repr(δx) norm_step=norm(δx, 2) α αmax repr(x) y fcalls norm_gprev=norm(∇f!(similar(x), x - α*δx)) norm_gnew=norm(∇) prod(diag(hess.U))

        # check for convergence
        convcond(x, xpre, y, ypre, ∇) && return NewtonResult(true, x, i, totfcalls)
    end
    return NewtonResult(false, x, maxiter, totfcalls)
end
