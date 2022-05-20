# CubicEoS.jl

[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)

CubicEoS.jl is an extensible package for thermodynamic calculations. The package defines an abstract interface needed to implement an equation of state (EoS) and uses it to calculate NVT phase equilibrium.

As an example of EoS, the package implements the general cubic equation of state [[Brusilovsky, SPE Reservoir Engineering, February 1992](https://doi.org/10.2118/20180-PA)]. To implement an EoS you need, use `src/BrusilovskyEoS` submodule as a template.

TODO: add links to other extensions: CPPCSAFT and MBWREoS.

What CubicEoS.jl do provide?

- Interface needed to solve NVT phase equilibrium problem;
- NVT solvers: phase stability and phase split which
  - Support automatic differentation;
  - Based on optimization: currently uses Cholesky BFGS with control of step from **[Downhill.jl](https://github.com/vvpisarev/Downhill.jl)**;
  - Provide choice of internal variables for better scaling of a problem;
- Implementation of the general equation of state.

Most of implemented functions

- are designed in a zero-allocating way, such functions has optional `buf` keyword;
- has a desctructive option `func!`.

## Installation

The package registry is **[LammpsToolsRegistry](https://github.com/stepanzh-test-org/LammpsToolsRegistry)**. So, you need to add the registry and then install CubicEoS.jl in a usual Julia-way.

```julia
]pkg> registry add https://github.com/stepanzh-test-org/LammpsToolsRegistry
pkg> add CubicEoS
```

## Documentation

Currently, only docstrings are provided. You may also take a look at tests with minimal working examples.

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
compressibility(mixture, nmol, pressure, RT, phase='g')
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
issuccess, isstable, results = vt_stability(mix, nmol, volume, RT, StateVariables)
```

**Phase split.** To calculate NVT phase equilibrium use

```julia
result = vt_split(mix, nmol, volume, RT, trial_concentration, StateVariables)
```

where type `StateVariables` defines an internal variables used at phase split stage in optimization solver (e.g. `CubicEoS.PhysicalState`). For the split `trial_concentration` is given from results of the stability.

There are several options to control behavior of the solvers. Check docstrings of `vt_stability` and `vt_split` for overview.

## Extensions to other equations of state

Currently, to implement an EoS you need, use `src/BrusilovskyEoS` submodule as a template (check functions that add methods to the main module, e.g. `CubicEoS.pressure`).

Briefly speaking, the interface can be divided into servicing and physical categories. The servicing interface mostly requires getter-like functions. The physical interface requires minimal set of functions to compute pressure, activity coefficients and constraints on phases introducing by an EoS.

Currently, there are several EoS under development:

TODO: add links to other extensions: CPPCSAFT and MBWREoS.

