module CubicEoS

export AbstractEoSComponent, AbstractEoSMixture
export molar_mass, carbon_number, name, describe, components, thermo_buffer, ncomponents
export pressure, wilson_saturation_pressure, compressibility
export log_c_activity!, log_c_activity, log_c_activity_wj!, log_c_activity_wj
export vt_stability, vt_stability!
export vt_split, vt_split!
export converged
export BrusilovskyEoSComponent, BrusilovskyEoSMixture

using CubicEoSDatabase
using Downhill
using LinearAlgebra
using StructArrays

include("constants.jl")
include("types.jl")
include("error.jl")
include("interface.jl")
include("solvecubic.jl")
include("basic_thermo.jl")
include("chempotential.jl")
include("vt_stability/vt_stability.jl")
include("vt_flash/vt_flash.jl")

include("BrusilovskyEoS/BrusilovskyEoS.jl")
using .BrusilovskyEoS

end # module
