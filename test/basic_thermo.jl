@testset "Basic thermo" begin
    C₁C₅ = CubicEoS.load(BrusilovskyEoSMixture; names=("methane", "n-pentane"))

    @testset "eos_parameters" begin
        nmol = [0.8, 1-0.8]
        RT = 300 * CubicEoS.GAS_CONSTANT_SI
        @testset "Calls w/wo buffers" begin
            buffer_st = thermo_buffer(C₁C₅)
            buffer_nt = (ai=similar(nmol), aij=zeros((length(nmol), length(nmol))))
            buffer_dict = Dict(:ai => similar(nmol), :aij => zeros((length(nmol), length(nmol))))
            # no buffer vs buffer struct
            @test CubicEoS.eos_parameters(C₁C₅, nmol, RT) == CubicEoS.eos_parameters(C₁C₅, nmol, RT; buf=buffer_st)
            # no buffer vs buffer namedtuple
            @test CubicEoS.eos_parameters(C₁C₅, nmol, RT) == CubicEoS.eos_parameters(C₁C₅, nmol, RT; buf=buffer_nt)
            # no buffer vs buffer dict
            @test CubicEoS.eos_parameters(C₁C₅, nmol, RT) == CubicEoS.eos_parameters(C₁C₅, nmol, RT; buf=buffer_dict)
        end
        @testset "Homogeneity" begin
            @test_skip CubicEoS.eos_parameters(C₁C₅, nmol, RT) == CubicEoS.eos_parameters(C₁C₅, 5 * nmol, RT)
        end
    end
end
