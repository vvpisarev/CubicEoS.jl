@testset "basic_thermo.jl" begin
    C₁C₅ = CubicEoS.load(BrusilovskyEoSMixture; names=("methane", "n-pentane"))
    RT = 300 * CubicEoS.GAS_CONSTANT_SI
    @testset "Calls w/wo buffers" begin
        buffer = thermo_buffer(C₁C₅)

        N = 2 * [0.8, 1-0.8]
        V = 2.5
        @test pressure(C₁C₅, N, V, RT) == pressure(C₁C₅, N, V, RT; buf=buffer)
    end
    @testset "Pressure" begin
        @testset "Homogeneity" begin
            V = 1
            χ = [0.8, 1-0.8]
            for ∑N in 0.1:0.2:1.9
                N = χ * ∑N
                υ = V / ∑N
                @test pressure(C₁C₅, N, V, RT) ≈ pressure(C₁C₅, χ, υ, RT)
            end
        end
    end
end
