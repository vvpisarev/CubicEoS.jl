@testset "chempotential.jl" begin
    C₁C₅ = CubicEoS.load(BrusilovskyEoSMixture; names=("methane", "n-pentane"))
    N = 1e-2 .+ rand(2)
    V = 1e-2 + rand()
    RT = 300 * CubicEoS.GAS_CONSTANT_SI

    @testset "func vs func! versions" begin
        vec = similar(N)
        log_c_activity!(vec, C₁C₅, N, V, RT)
        @test vec == log_c_activity(C₁C₅, N, V, RT)

        matr = vec * vec'
        log_c_activity_wj!(vec, matr, C₁C₅, N, V, RT)
        @test (vec, matr) == log_c_activity_wj(C₁C₅, N, V, RT)
    end

    @testset "Buffered versions" begin
        buf = thermo_buffer(C₁C₅)
        @test log_c_activity(C₁C₅, N, V, RT) == log_c_activity(C₁C₅, N, V, RT; buf=buf)
        @test log_c_activity_wj(C₁C₅, N, V, RT) == log_c_activity_wj(C₁C₅, N, V, RT; buf=buf)
    end

    @testset "Homogeneity" begin
        for α in 0.1:0.2:1.9
            @test log_c_activity(C₁C₅, N, V, RT) ≈ log_c_activity(C₁C₅, α .* N, α * V, RT)

            v, m = log_c_activity_wj(C₁C₅, N, V, RT)
            vα, mα = log_c_activity_wj(C₁C₅, α .* N, α * V, RT)
            @test v ≈ vα
            # m is jacobian
            # while activity has 0-order homogeneity,
            # its derivatives (jacobian) is -1-order homogeneous
            @test m ≈ α .* mα
        end
    end
end
