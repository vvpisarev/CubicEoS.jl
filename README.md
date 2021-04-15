# CubicEoS.jl
The package implements functions to work with substances and mixtures described by cubic equations of state.

So far, the general cubic equation of state [[Brusilovsky, SPE Reservoir Engineering, February 1992](https://doi.org/10.2118/20180-PA)] is implemented.

## Brusilovsky equation of state

To construct an set of parameters for a component, use
```julia
BrusilovskyEoSComponent(;name="No Name", critical_pressure, critical_temperature, acentric_factor, Omegac, Zc, Psi, molar_mass, carbon_number::Integer)
```
The values of parameters `Omegac`, `Zc` and `Psi` for hydrocarbons may be found in Brusilovsky's paper (https://doi.org/10.2118/20180-PA).

The temperatures must be in absolute scale (e.g., in Kelvins).

The parameters can be loaded from a file (using **CubicEoSDatabase.jl**):
```julia
methane = load(BrusilovskyEoSComponent, name = "methane", physics_db = "martinez.csv", eos_db = "brusilovsky.csv")
```

Mixtures are constructed via
```julia
BrusilovskyEoSMixture(; components::AbstractVector{<:BrusilovskyEoSComponent}, constant, linear, quadratic)
```
where `constant`, `linear` and `quadratic` are matrices of constant, linear and quadratic in temperature terms for Zudkevitch-Joffe corrections.

The parameters can be loaded from a file (using **CubicEoSDatabase.jl**):
```julia
c1c5 = load(BrusilovskyEoSMixture, names = ["methane", "n-pentane", comp_physics_db = "martinez.csv", comp_eos_db = "brusilovsky.csv", mix_eos_db = "brusilovsky_mix.csv")
```

## Basic thermodynamics

To get the pressure of a pure component:
```julia
pressure(component; nmol, volume, temperature)
```

To get the estimate of the saturation pressure at a given temperature:
```julia
wilson_saturation_pressure(component, RT)
```

To get compressibility at given pressure and temperature:
```julia
compressibility(mixture, molar_composition, pressure, RT, phase = 'g'
```

## Phase stability

To check if a single-phase state is stable:
```julia
vt_stability(mix::BrusilovskyEoSMixture, nmol::AbstractVector, volume, RT)
```
