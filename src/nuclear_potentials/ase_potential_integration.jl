# Integrations with ASEPotential.jl

struct ASEPotentialParameters <: NuclearPotentialParameters
    element::ElementCoulomb  # Why? Could be a generic element (it's actually just by accident DFTK returns an ElementCoulomb here. This is only because we don't know better, but this may change in the future)
    lattice::AbstractArray{Quantity, 2}
    parameters::ASEPotential.ASECalculatorParameters
end

function forces(state::AtomicConfiguration, parameters::ASEPotentialParameters)
    atoms = ASEAtoms(state, parameters.element, parameters.lattice).atoms
    forces = get_forces(atoms, parameters.parameters)
    SVector{3}.(eachrow(forces))
end

function potential_energy(state::AtomicConfiguration, parameters::ASEPotentialParameters)
    atoms = ASEAtoms(state, parameters.element, parameters.lattice).atoms
    get_potential_energy(atoms, parameters.parameters)
end
