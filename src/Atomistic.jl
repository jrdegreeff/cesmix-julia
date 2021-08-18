module Atomistic

using ASEPotential
using DFTK
using NBodySimulator
using Plots
using PyCall
using StaticArrays
using Unitful
using UnitfulAtomic
using UnitfulRecipes


include("exceptions.jl")

export MassBodies, DFTKAtoms, ASEAtoms
include("bodies.jl")

export NuclearPotentialParameters, forces, potential_energy, LJParameters
include("nuclear_potentials/abstractions.jl")
export ASEPotentialParameters, ASEPotentialParameters
include("nuclear_potentials/ase_potential_integration.jl")
export DFTKParameters, dftk_atoms, analyze_convergence
include("nuclear_potentials/dftk_integration.jl")

# Personally I'm not a fan of long line lengths, but everyone is different.
export MolecularDynamicsParameters, MolecularDynamicsResult, simulate, get_bodies, get_time_range, plot_temperature, plot_temperature!, plot_energy, plot_energy!, plot_rdf, calculate_rdf
include("molecular_dynamics/abstractions.jl")
export NBSParameters, NBSResult
include("molecular_dynamics/nbs_integration.jl")

export write_nbs_animation, write_ase_trajectory
include("io/visualization.jl")

end
