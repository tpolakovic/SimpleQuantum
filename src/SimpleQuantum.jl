module SimpleQuantum

using LinearAlgebra
using Unitful
import Unitful: Å, eV
using RangeHelpers: range, around
using Match
using Colors
using GLMakie
using Accessors

include("crystal.jl")
include("interface.jl")
include("tb.jl")
include("nfe.jl")
include("misc.jl")

# Unit definitons

a₀ = 1.889726125Å
Ha = 27.2eV

export a₀,
    Ha,
    Å,
    eV,
    Lattice,
    UnitCell,
    Crystal,
    plotcrystal!,
    fermilevel,
    shiftenergy,
    evals,
    evecs,
    kvecs,
    solve,
    ReciprocalPath,
    plotSolution,
    DOS,
    Hoppings,
    addhop!,
    addonsite!,
    addoverlap!,
    TightBindingHamiltonian,
    PseudoPotentialHamiltonian
end
