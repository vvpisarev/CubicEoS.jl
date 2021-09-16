@testset "solvecubic.jl" begin
    mix = load(BrusilovskyEoSMixture; names=("methane",))
    solve_cubic = CubicEoS.solve_cubic

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
end
