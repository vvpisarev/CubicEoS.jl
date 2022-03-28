struct ConstrainStepZeroDirectionError <: Exception
    index::Int
    statevalue::Float64
end

Base.showerror(io::IO, e::ConstrainStepZeroDirectionError) = begin
    println(io, "ConstrainStepZeroDirectionError along direction ", e.index, ":")
    print(io, "illegal state value ", e.statevalue)
end

struct ConstrainStepLowerBoundError{V1,V2} <: Exception
    x::V1
    dir::V2
end

Base.showerror(io::IO, e::ConstrainStepLowerBoundError) = begin
    println(io, "ConstrainStepLowerBoundError")
    println(io, "x = ", repr(e.x))
    print(io, "direction = ", repr(e.dir))
end
