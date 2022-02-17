# CubicEoS.jl

[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)

CubicEoS.jl implements functions to work with substances and mixtures described by cubic equations of state. This includes

- Basic thermodynamics: pressure, z-factor, chemical potential etc.;
- Check of NVT stability of single-phase state;
- Phase equilibrium calculation for NVT variables.

So far, the general cubic equation of state [[Brusilovsky, SPE Reservoir Engineering, February 1992](https://doi.org/10.2118/20180-PA)] is implemented.

Most of implemented functions

- are designed in a zero-allocating way, such functions has optional `buf` keyword;
- has a desctructive option `func!`.

## Installation

The package registry is **[LammpsToolsRegistry](https://github.com/stepanzh-test-org/LammpsToolsRegistry)**. So, you need to add the registry and then install CubicEoS.jl in a usual Julia-way.

```julia
]pkg> registry add https://github.com/stepanzh-test-org/LammpsToolsRegistry
pkg> add CubicEoS
```

## Brusilovsky equation of state

Components and mixtures may be constructed either explicitly or by loading.

A **component** can be defined explicitly
```julia
BrusilovskyEoSComponent(; name="No Name", critical_pressure, critical_temperature, acentric_factor, Omegac, Zc, Psi, molar_mass, carbon_number::Integer)
```
where the values of parameters `Omegac`, `Zc` and `Psi` may be found in Brusilovsky's paper (https://doi.org/10.2118/20180-PA).

The temperatures must be in absolute scale (e.g., in Kelvins).

Or, the component can be loaded from a file (using **[CubicEoSDatabase.jl](https://github.com/stepanzh/CubicEoSDatabase.jl)**)
```julia
methane = CubicEoS.load(BrusilovskyEoSComponent; name="methane"[, custom_databases...])
```

A **mixture** is constructed via

```julia
BrusilovskyEoSMixture(; components::AbstractVector{<:BrusilovskyEoSComponent}, constant, linear, quadratic)
```
where `constant`, `linear` and `quadratic` are matrices of constant, linear and quadratic in temperature terms for Zudkevitch-Joffe corrections.

Or can be loaded from a file (using **[CubicEoSDatabase.jl](https://github.com/stepanzh/CubicEoSDatabase.jl)**):
```julia
c1c5 = CubicEoS.load(BrusilovskyEoSMixture; names=("methane", "n-pentane")[, custom_databases...])
```

## Basic thermodynamics

Basic thermodynamics includes pressure, Wilson saturation pressure and z-factor (`compressibility`).

```julia
pressure(component or mixture, nmol, volume, temperature)
```

```julia
wilson_saturation_pressure(component, RT)
```

```julia
compressibility(mixture, molar_composition, pressure, RT, phase='g')
```

## Chemical potential

The packages includes functions for calculating activity coefficient and its Jacobian matrix for a mixture defined by Brusilovsky EoS.

```julia
log_ca = log_c_activity(mixture, nmol, volume, RT)
log_ca, jac = log_c_activity_wj(mixture, nmol, volume, RT)
```

In case of a component you may use a mixture of one component.

## NVT phase equilibrium

**Phase stability.** To check if a single-phase state is stable, defined in NVT variables, use

```julia
isstable, vtstab_results = vt_stability(mix, nmol, volume, RT)
```

**Flash.** To calculate NVT phase equilibrium use

```julia
# quasi-Newton phase split
flash_result = vt_flash(mix, nmol, volume, RT, StateVariables[; tol, maxiter])
# Newton phase split
flash_result = vt_flash_newton(mix, nmol, volume, RT, StateVariables[; tol, maxiter])
```

where type `StateVariables` defines an internal variables used at phase split stage in optimization solver (e.g. `CubicEoS.PhysicalState`). Quasi-Newton solver is implemented in **[Downhill.jl](https://github.com/vvpisarev/Downhill.jl)**.

