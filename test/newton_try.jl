using CubicEoS
using LinearAlgebra


function quadratic(M::AbstractMatrix)
    @assert isposdef(M)
    @assert abs(det(M)) > 0

    F(x) = x' * M * x
    dF!(g, x) = begin g .= (M' + M) * x; g end
    ddF!(H, x) = begin H .= (M' + M); H end

    return F, dF!, ddF!
end

function quadratic(M::AbstractMatrix, b::AbstractVector)
    @assert isposdef(M)
    @assert abs(det(M)) > 0

    F(x) = x' * M * x - 2 * x' * b
    dF!(g, x) = begin g .= (M' + M) * x .- 2 .* b; g end
    ddF!(H, x) = begin H .= (M' + M); H end

    return F, dF!, ddF!
end

function system_test(M::AbstractMatrix, b::AbstractVector)
    f, df, ddf = quadratic(M, b)
    x₀ = 1000 .* rand(size(M, 1))
    optresult = CubicEoS.newton(f, df, ddf, x₀;
        maxiter=100,
        gtol=1e-3,
        constrain_step=(x, δ)->1.0,
    )
    dump(optresult)
    @show norm(optresult.argument, 2)

    @show f(optresult.argument)

    @show M \ b
    @show f(M \ b)
end

function quadratic_test(M::AbstractMatrix)
    f, df, ddf = quadratic(M)

    x₀ = 1000 .* rand(size(M, 1))
    optresult = CubicEoS.newton(f, df, ddf, x₀;
        maxiter=100,
        gtol=1e-6,
        constrain_step=(x, δ)->1.0,
    )
    dump(optresult)
    @show norm(optresult.argument, 2)
end

function quadratic_test(n::Integer=20)
    A = rand(n, n)
    return quadratic_test(A*A')
end

# quadratic_test(40)
# let A = [1 2; 4 5]
#     system_test(A*A', [1, 1])
# end
let n = 10, A = rand(n, n)
    system_test(A*A', rand(n))
end

