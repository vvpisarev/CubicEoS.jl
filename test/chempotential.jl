@testset "chempotential.jl" begin
    C₁C₅ = CubicEoS.load(BrusilovskyEoSMixture; names=("methane", "n-pentane"))
    N = 1e-2 .+ rand(2)
    V = 1e-2 + rand()
    RT = 300 * CubicEoS.GAS_CONSTANT_SI


    @testset "func vs func! versions" begin
        vec = similar(N)
        log_c_activity!(vec, C₁C₅, N, V, RT)
        @test vec == log_c_activity(C₁C₅, N, V, RT)
    end

    @testset "Buffered versions" begin
        buf = thermo_buffer(C₁C₅)
        @test log_c_activity(C₁C₅, N, V, RT) == log_c_activity(C₁C₅, N, V, RT; buf = buf)
    end

    @testset "Homogeneity" begin
        for α in 0.1:0.2:1.9
            @test log_c_activity(C₁C₅, N, V, RT) ≈ log_c_activity(C₁C₅, α .* N, α * V, RT)
        end
    end
end
