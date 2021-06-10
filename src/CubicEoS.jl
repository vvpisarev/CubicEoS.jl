module CubicEoS

export AbstractEoSComponent, AbstractEoSMixture
export BrusilovskyEoSComponent, BrusilovskyEoSMixture, thermo_buffer
export molar_mass, carbon_number, name, describe, load, components # Do we need to export `load`?
export ncomponents, pressure, wilson_saturation_pressure, compressibility
export log_c_activity!, log_c_activity, log_c_activity_wj!, log_c_activity_wj
export vt_stability, vt_stability_buffer

include("constants.jl")
include("types.jl")
include("interface.jl")
include("dbload.jl")
include("basic_thermo.jl")
include("chempotential.jl")
include("vt_stability.jl")

end # module
