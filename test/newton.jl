@testset "newton.jl" begin
    booth(x) = (x[1] + 2*x[2] - 7)^2 + (2*x[1] + x[2] - 5)^2
    gradbooth!(g, x) = begin g .= Float64[10*x[1] + 8*x[2] - 34, 8*x[1] + 10*x[2] - 38]; g end
    hessbooth!(hess, x) = begin hess .= Float64[10 8; 8 10]; hess end

    xoptim = [1, 3]  # booth function minimum
    x = [4, 2]

    @testset "Stop by gradient" begin
        ptrx = pointer(x)
        xbefore = copy(x)
        result = CubicEoS.newton(booth, gradbooth!, hessbooth!, x)
        @testset "Sameness of argument" begin
            @test pointer(x) == ptrx
            @test x == xbefore
        end

        @testset "Converged test" begin
            @test result.converged
            @test result.argument ≈ xoptim
            @test booth(result.argument) ≈ 0 atol=1e-10
        end
    end

    @testset "Don't stop" begin
        result = CubicEoS.newton(booth, gradbooth!, hessbooth!, x; gtol=NaN, maxiter=2)
        @test !result.converged
    end

    @testset "Custom stop condition" begin
        result = CubicEoS.newton(booth, gradbooth!, hessbooth!, x;
                    gtol=NaN,
                    maxiter=2,
                    convcond=(x, xpre, y, ypre, g) -> abs(y - ypre) ≤ 1e-6,
        )
        @test result.converged
        @test result.argument ≈ xoptim
    end
end
