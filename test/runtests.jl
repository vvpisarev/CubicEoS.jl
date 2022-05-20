using CubicEoS
using CubicEoSDatabase: Data
using Test

"Non-recursively compare field values of `x` and `y`."
function compare_structs(x::T, y::T) where {T}
    for field in fieldnames(T)
        @test getfield(x, field) == getfield(y, field)
    end
end

@testset "CubicEoS.jl" begin
    include("solvecubic.jl")
    include("newton.jl")
    include("BrusilovskyEoS/BrusilovskyEoS.jl")
    include("vt_stability/vt_stability.jl")
    include("vt_split/vt_split.jl")
end
