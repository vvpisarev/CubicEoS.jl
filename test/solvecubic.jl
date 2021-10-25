@testset "solvecubic.jl" begin
    mix = load(BrusilovskyEoSMixture; names=("methane",))
    solve_cubic = CubicEoS.solve_cubic

    @testset "Unstable 1-phase state" begin
        # definitely inside "bulk" of two-phase region
        @test !first(vt_stability(mix, [0.01], 1e-6, CubicEoS.GAS_CONSTANT_SI * 140.0))
    end

    @testset "Degenerate equation" begin
        @test first(vt_stability(mix, [0.0109], 1e-6, 1588.0623600672689))
        @test first(vt_stability(mix, [0.0109], 1e-6, CubicEoS.GAS_CONSTANT_SI * 190.6))
    end

    @testset "Triple root" begin
        @test solve_cubic(1, -3, 3, -1) == (1.0, 1.0, 1.0)
        @test solve_cubic(1, 3, 3, 1) == (-1.0, -1.0, -1.0)
    end

    @testset "One real root" begin
        @test let
            p, q, r = solve_cubic(1, -5, 1, -5)
            p ≈ 5.0 && isnan(q) && isnan(r)
        end
        @test let
            p, q, r = solve_cubic(1, -5, 9, -5)
            p ≈ 1.0 && isnan(q) && isnan(r)
        end
    end

    @testset "solve_cubic!" begin
        solve_cubic! = CubicEoS.solve_cubic!
        @test let
            coefs = 1, -5, 1, -5
            p, q, r = solve_cubic(coefs...)
            buf = Vector{Float64}(undef, 5)
            solve_cubic!(buf, coefs...)
            any(buf[1:3] .== (p, q, r))
        end
    end
end
