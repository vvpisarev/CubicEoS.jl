include("types.jl")
include("nvt.jl")
# include("state_physical.jl")
# include("state_ratio.jl")
# include("state_idealidentity.jl")
# include("newton.jl")

# """
#     vt_flash(mix, nmol, volume, RT, StateVariables[; tol, chemtol=tol, presstol=tol, maxiter=100])

# Two-phase thermodynamical equilibrium solver for `mix`ture at given moles `nmol`, `volume`
# and thermal energy `RT` (VT-flash).
# Includes two stages, the first is stability checking of single-phase state,
# if the state is unstable, then an initial two-phase state is constructed,
# and phase-split is performed with `StateVariables` using Cholesky-BFGS optimization.

# For two-phase state the equilibrium is considered, when

# 1. Chemical potentials are equal in a sense

# ```
#  1
# --- maxᵢ |μᵢ' - μᵢ''| < chemtol
#  RT
# ```

# 2. Pressures are equals in a sense
# ```
# |P' - P''| volume
# ----------------- < presstol,
#    RT sum(nmol)
# ```
# where `i` is component index, and `'`, `''` are indexes of phases.

# Return [`VTFlashResult`](@ref).

# See also [`CubicEoS.vt_flash!`](@ref), [`vt_flash_newton`](@ref).

# # Arguments

# - `mix::BrusilovskyEoSMixture{T}`: mixture;
# - `nmol::AbstractVector`: moles of mixtures' components [mole];
# - `volume::Real`: volume of mixture [meter³];
# - `RT::Real`: use to specify temperature of mixture,
#     `CubicEoS.GAS_CONSTANT_SI * temperature`, [Joule / mole];
# - `StateVariables::Type{<:AbstractVTFlashState}`: one of state variables to use internally
#     in phase-split stage.

# # Optional arguments

# - `chemtol::Real=tol`: tolerance for chemical potentials of components;
# - `presstol::Real=tol`: tolerance for pressures of phases;
# - `tol::Real=1024*eps(T)`: default tolerance for both `chemtol` and `presstol`,
#     `T` is AbstractFloat type defined by `mixture`'s type;
# - `maxiter::Integer`: maximum allowed steps in phase-split stage.
# """
# function vt_flash(
#     mix::BrusilovskyEoSMixture{T},
#     nmol::AbstractVector,
#     volume::Real,
#     RT::Real,
#     StateVariables::Type{<:AbstractVTFlashState};
#     tol::Real=1024*eps(T),
#     chemtol::Real=tol,
#     presstol::Real=tol,
#     maxiter::Integer=100,
# ) where {T}
#     stabconverged, singlephase, stability_tries = vt_stability(mix, nmol, volume, RT)

#     !stabconverged && error("VT-Stability does not converged")

#     if singlephase
#         concentration = nmol ./ volume
#         saturation = 1
#         state = StateVariables(concentration, saturation, nmol, volume)
#         return __vt_flash_single_phase_result(state, mix, nmol, volume, RT)
#     end

#     concentration = __vt_flash_init_conc_choose(stability_tries)
#     saturation = __find_saturation_negative_helmdiff(mix, nmol, volume, RT, concentration;
#         maxsaturation=0.25,
#         maxiter=50,
#         scale=0.5,
#         helmdifftol=-1e-7/RT,
#     )

#     if isnan(saturation)
#         @error "VTFlash: Initial state was not found!" mixture=mix nmol=repr(nmol) volume=volume RT=RT
#         error("VTFlash: Initial state was not found!")
#     end

#     state = StateVariables(concentration, saturation, nmol, volume)

#     return vt_flash!(state, mix, nmol, volume, RT;
#         chemtol=chemtol,
#         presstol=presstol,
#         maxiter=maxiter
#     )
# end

# """
#     vt_flash!(unstable_state, mix, nmol, volume, RT; chemtol, presstol, maxiter)

# Perform split phase of VT-flash from an `unstable_state::AbstractVTFlashState`,
# which will be destructed.

# For rest of arguments see [`vt_flash`](@ref).

# Return [`VTFlashResult`](@ref).
# """
# function vt_flash!(
#     unstable_state::AbstractVTFlashState,
#     mix::BrusilovskyEoSMixture,
#     nmol::AbstractVector,
#     volume::Real,
#     RT::Real;
#     chemtol::Real,
#     presstol::Real,
#     maxiter::Int,
# )
#     state = unstable_state
#     statex = value(state)

#     # initial hessian
#     hessian = Matrix{Float64}(undef, (size(statex, 1), size(statex, 1)))
#     hessian = hessian!(hessian, state, mix, nmol, volume, RT)

#     constrain_step, helmgrad!, helmdiff! = __vt_flash_optim_closures(
#         state, mix, nmol, volume, RT
#     )

#     convcond = __convergence_closure(state, mix, nmol, volume, RT;
#         chemtol=chemtol,
#         presstol=presstol,
#     )

#     # run optimize
#     optmethod = Downhill.CholBFGS(statex)
#     Downhill.reset!(optmethod, statex, hessian)
#     optimresult = Downhill.optimize!(helmdiff!, optmethod, statex;
#         gtol=NaN,
#         convcond=convcond,
#         maxiter=maxiter,
#         constrain_step=constrain_step,
#         reset=false,
#     )
#     statex .= optimresult.argument

#     return __vt_flash_two_phase_result(state, mix, nmol, volume, RT, optimresult)
# end

# function __vt_flash_single_phase_result(
#     state::S,
#     mix::BrusilovskyEoSMixture,
#     nmol::AbstractVector{T},
#     volume::Real,
#     RT::Real
# ) where {S, T}
#     return VTFlashResult{T, S}(;
#         converged=true,
#         singlephase=true,
#         RT=RT,
#         nmolgas=nmol,
#         volumegas=volume,
#         nmolliq=similar(nmol),
#         volumeliq=0,
#         state=state,
#         iters=-1,
#         calls=-1,
#     )
# end

# function __find_saturation_negative_helmdiff(
#     mix::BrusilovskyEoSMixture,
#     nmolb::AbstractVector,
#     volumeb::Real,
#     RT::Real,
#     concentration::AbstractVector;
#     maxsaturation::Real,
#     maxiter::Int,
#     scale::Real,
#     helmdifftol::Real,
#     thermo_buf::BrusilovskyThermoBuffer=thermo_buffer(mix),
# )
#     # Δa = a' + a'' - abase
#     abase = helmholtz(mix, nmolb, volumeb, RT; buf=thermo_buf)

#     nmol1 = similar(nmolb, Float64)
#     nmol2 = similar(nmolb, Float64)

#     function helmdiff(saturation::Real)
#         volume1 = volumeb * saturation
#         @. nmol1 = volume1 * concentration
#         a1 = helmholtz(mix, nmol1, volume1, RT; buf=thermo_buf)

#         volume2 = volumeb - volume1
#         @. nmol2 = nmolb - nmol1
#         a2 = helmholtz(mix, nmol2, volume2, RT; buf=thermo_buf)
#         return (a1 + a2) - abase
#     end

#     saturation = float(maxsaturation)
#     for i in 1:maxiter
#         try
#             Δa = helmdiff(saturation)
#             Δa < - abs(helmdifftol) && return saturation
#         catch e
#             isa(e, UndefVarError) && throw(e)  # syntax
#         end
#         saturation *= scale
#     end
#     return NaN
# end

# function __convergence_closure(
#     state1::AbstractVTFlashState,
#     mix::BrusilovskyEoSMixture,
#     nmolb::AbstractVector,
#     volumeb::Real,
#     RT::Real;
#     chemtol::Real,
#     presstol::Real,
#     buf::BrusilovskyThermoBuffer=thermo_buffer(mix),
# )
#     state1x = value(state1)

#     g1 = Vector{Float64}(undef, ncomponents(mix) + 1)
#     g2 = similar(g1)
#     diff = similar(g1)
#     nmol2 = similar(g1, Float64, ncomponents(mix))

#     function convcond(x::V, xpre::V, y::T, ypre::T, g::V) where {T<:Real,V<:AbstractVector}
#         state1x .= x
#         nmol1, vol1 = nmolvol(state1, nmolb, volumeb)
#         g1 = CubicEoS.nvtgradient!(g1, mix, nmol1, vol1, RT; buf=buf)

#         @. nmol2 = nmolb - nmol1
#         vol2 = volumeb - vol1
#         g2 = CubicEoS.nvtgradient!(g2, mix, nmol2, vol2, RT; buf=buf)

#         #=
#                    [μ₁' - μ₁'', ..., μₙ' - μₙ''; -P' + P'']ᵀ
#             diff = ----------------------------------------
#                                     RT
#         =#
#         @. diff = g1 - g2

#         condchem = let diffchem = (@view diff[1:end-1])
#             # max_i |chem'_i - chem''_i| / RT < tolerance
#             maximum(abs, diffchem) < chemtol
#         end
#         condpress = let diffpress = diff[end]
#             # |P' - P''| V
#             # ------------ < tolerance
#             #    RT ΣNᵢ
#             abs(diffpress) * volumeb / sum(nmolb) < presstol
#         end

#         return condchem && condpress
#     end

#     return convcond
# end

# function __sort_phases!(mix, nmol₁, V₁, nmol₂, V₂, RT)
#     P₁ = pressure(mix, nmol₁, V₁, RT)  # they should be equal
#     P₂ = pressure(mix, nmol₂, V₂, RT)

#     Z₁ = P₁ * V₁ / (sum(nmol₁) * RT)  # seems can be reduced to Vᵢ / sum(nmolᵢ)
#     Z₂ = P₂ * V₂ / (sum(nmol₂) * RT)

#     if Z₂ > Z₁  # □₂ is gas state, need exchange
#         P₁, P₂ = P₂, P₁
#         Z₁, Z₂ = Z₂, Z₁
#         V₁, V₂ = V₂, V₁

#         for i in eachindex(nmol₁, nmol₂)
#             nmol₁[i], nmol₂[i] = nmol₂[i], nmol₁[i]
#         end
#     end
#     # now □₁ is gas, □₂ is liquid
#     return nmol₁, V₁, nmol₂, V₂
# end

# """
# Extracts vt-state from `optresult` (Downhill obj).
# Sorts variables into gas and liquid.
# Returns corresponding `VTFlashResult`.
# """
# function __vt_flash_two_phase_result(
#     state::S,
#     mix::BrusilovskyEoSMixture{T},
#     nmol::AbstractVector{T},
#     volume::Real,
#     RT::Real,
#     optresult,
# ) where {T, S<:AbstractVTFlashState}
#     nmol1, volume1 = nmolvol(state, nmol, volume)

#     nmol2 = nmol .- nmol1
#     volume2 = volume - volume1

#     nmolgas, volgas, nmolliq, volliq = __sort_phases!(mix, nmol1, volume1, nmol2, volume2, RT)

#     return VTFlashResult{T, S}(;
#             converged=optresult.converged,
#             singlephase=false,
#             RT=RT,
#             nmolgas=nmolgas,
#             volumegas=volgas,
#             nmolliq=nmolliq,
#             volumeliq=volliq,
#             state=state,
#             iters=optresult.iterations,
#             calls=optresult.calls,
#     )
# end


# "Return concentration of state with minimum energy from vt-stability tries."
# function __vt_flash_init_conc_choose(vt_stab_tries)
#     Dmin = Inf
#     index_min = -1
#     for (i, state) in enumerate(vt_stab_tries)
#         if state.issuccess && !state.isstable && state.energy_density < Dmin
#             index_min = i
#             Dmin = state.energy_density
#         end
#     end
#     index_min == -1 && error("Stability tries are inconsistent: Can't choose the one with the lowest energy")
#     return vt_stab_tries[index_min].concentration
# end
