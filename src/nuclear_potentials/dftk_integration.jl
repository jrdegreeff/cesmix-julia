# Integrations with DFTK.jl

Base.@kwdef struct DFTKParameters <: NuclearPotentialParameters
    # Here I would only take a basis and a kwargs object passed to the SCF. That's it.
    psp::ElementPsp
    lattice::AbstractArray{Quantity, 2}
    Ecut::Quantity   # This is a great example why putting too many types is bad. DFTK is quite flexible wrt. Ecut (integers, floats, unitful, all allowed). By specifying ::Quantity here you *reduce* the flexibility of your interface.
    kgrid::AbstractVector{Integer}
    n_bands::Union{Integer, Nothing} = nothing
    tol::Union{AbstractFloat, Nothing} = nothing
    damping::Union{AbstractFloat, Nothing} = nothing
    mixing::Union{DFTK.Mixing, Nothing} = nothing
    previous_scfres::Base.RefValue{Any} = Ref{Any}()
end
DFTKParameters(parameters::DFTKParameters, Ecut::Quantity) = DFTKParameters(
    psp=parameters.psp,
    lattice=parameters.lattice,
    Ecut=Ecut,
    kgrid=parameters.kgrid,
    n_bands=parameters.n_bands,
    tol=parameters.tol,
    damping=parameters.damping,
    mixing=parameters.mixing,
)

function forces(state::AtomicConfiguration, parameters::DFTKParameters)
    # This is not the generic way you get forces in DFTK, this is specific for one atom only.
    compute_forces_cart(calculate_scf(state, parameters))[1]
end

function potential_energy(state::AtomicConfiguration, parameters::DFTKParameters)
    # You should check the SCF is converged. This is not guaranteed.
    calculate_scf(state, parameters).energies.total
end

function calculate_scf(state::AtomicConfiguration, parameters::DFTKParameters)
    setup_threading(n_blas=4)  # This should be done by the user.
    
    state = DFTKAtoms(state, parameters.psp, parameters.lattice)
    model = model_LDA(state.lattice, state.atoms)
    basis = PlaneWaveBasis(model, parameters.Ecut; kgrid=parameters.kgrid)

    extra_args = isassigned(parameters.previous_scfres) ? (ψ=parameters.previous_scfres[].ψ, ρ=parameters.previous_scfres[].ρ) : (; )
    scfres = self_consistent_field(basis; extra_args..., (f=>getfield(parameters, f) for f ∈ (:n_bands, :tol, :damping, :mixing) if getfield(parameters, f) !== nothing)...)
    parameters.previous_scfres[] = scfres
end

# This does not belong here. Actually this more belongs in DFTK as a utility. Nice idea to have this function, though!
function analyze_convergence(state::AtomicConfiguration, parameters::DFTKParameters, cutoffs::AbstractVector{<:Quantity})
    energies = Vector{Float64}()
    for Ecut ∈ cutoffs
        params = DFTKParameters(parameters, Ecut)
        scfres = calculate_scf(state, params)
        push!(energies, scfres.energies.total)
    end
    
    plot(
        title="DFTK Analysis",
        xlab="Ecut",
        ylab="Total Energy",
        legend=false,
        cutoffs,
        energies * u"hartree"
    )
end
