# Based on this guide: https://ase.tufts.edu/chemistry/lin/images/FortranMD_TeachersGuide.pdf
# Uses DFTK in place of LJ for the production stage as a proof-of-concept
# Note that the choice of parameters is for demonstration purposes and the results are non-physical

include("../../src/molecular_simulation.jl")
include("../../src/dftk_integration.jl")
include("../../src/nbs_extensions.jl")

include("./nbs_argon.jl")

N = 8
σ = auconvert(0.34u"nm")
box_size = 4σ # arbitrarly choosing 4σ

reference_temp = auconvert(94.4u"K")
thermostat_prob = 0.1 # this number was chosen arbitrarily

eq_steps = 20000
Δt = auconvert(1e-2u"ps")

eq_result, eq_bodies = argon_simulate_equilibration(N, box_size, Δt, eq_steps, reference_temp, thermostat_prob)

eq_stride = eq_steps ÷ 200

display(plot_temperature(eq_result, eq_stride))
display(plot_energy(eq_result, eq_stride))
display(plot_rdf(eq_result, sample_fraction=2))

# This scipt only runs DFTK once as a proof of concept -- note that this is not sensible for an ab initio simulation

dftk_parameters = DFTKForceGenerationParameters(
    box_size=box_size,
    psp=ElementPsp(:Ar, psp=load_psp(list_psp(:Ar, functional="lda")[1].identifier)),
    lattice=box_size * [[1. 0 0]; [0 1. 0]; [0 0 1.]],
    Ecut=10u"hartree",
    kgrid=[1, 1, 1],
    α=0.7,
    mixing=LdosMixing()
)

dftk_force_steps = 100

result, bodies = simulate(eq_bodies, dftk_parameters, box_size, Δt, dftk_force_steps)

dftk_force_stride = dftk_force_steps ÷ 10

# Ploting on separate plots because the timespan is so much smaller than in the first phase

display(plot_temperature(result, dftk_force_stride))
display(plot_energy(result, dftk_force_stride))
display(plot_rdf(result, σ=σ, sample_fraction=1))

;
