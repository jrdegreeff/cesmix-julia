# Nuclear Potential Abstractions

abstract type NuclearPotentialParameters end

# I would not provide these dummies, rather document what should be implemented.
# If not done, then julia ill throw.
#
# Interesting my feeling would have been parameters as first arg and state as second,
# but of course it does not really matter
function forces(state::AtomicConfiguration, parameters::NuclearPotentialParameters)
    throw(UnimplementedError(:forces, parameters))
end

function potential_energy(state::AtomicConfiguration, parameters::NuclearPotentialParameters)
    throw(UnimplementedError(:potential_energy, parameters))
end

Base.@kwdef struct LJParameters <: NuclearPotentialParameters
    ϵ::Quantity  # Again I would not put a single type here
    σ::Quantity
    R::Quantity
end
