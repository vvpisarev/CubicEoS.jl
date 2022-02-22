using CubicEoS
using CubicEoSDatabase
using Test


#= Utilities =#

"Non-recursively compare field values of `x` and `y`."
function compare_structs(x::T, y::T) where {T}
    for field in fieldnames(T)
        @test getfield(x, field) == getfield(y, field)
    end
end


@testset "CubicEoS.jl" begin
    include("dbload.jl")
    include("basic_thermo.jl")
    include("chempotential.jl")
    include("solvecubic.jl")
    include("vt_stability.jl")
    include("vt_flash/nvt.jl")
    include("vt_flash/state_idealidentity.jl")
    include("vt_flash/newton.jl")
    include("vt_flash/vt_flash.jl")
end
