@testset "Basic thermo" begin
    C₁C₅ = CubicEoS.load(BrusilovskyEoSMixture; names=("methane", "n-pentane"))
    RT = 300 * CubicEoS.GAS_CONSTANT_SI
    @testset "Calls w/wo buffers" begin
        nc = ncomponents(C₁C₅)
        buffer_st = thermo_buffer(C₁C₅)
        buffer_nt = (ai=zeros(nc), aij=zeros((nc, nc)))
        buffer_dict = Dict(:ai => zeros(nc), :aij => zeros((nc, nc)))
        # the following tests are grouped by
        # 1. no buffer vs buffer struct
        # 2. no buffer vs buffer namedtuple
        # 3. no buffer vs buffer dict
        N = 2 * [0.8, 1-0.8]
        V = 2.5
        for buffer in (buffer_st, buffer_nt, buffer_dict)
            @test CubicEoS.eos_parameters(C₁C₅, N, RT) == CubicEoS.eos_parameters(C₁C₅, N, RT; buf=buffer)
            @test pressure(C₁C₅, N, V, RT) == pressure(C₁C₅, N, V, RT; buf=buffer)
        end
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
