module SimpleQuantum

using LinearAlgebra
using Unitful
import Unitful: Å, eV
using RangeHelpers: range, around
using StaticArrays
using SplitApplyCombine
using Colors
using GLMakie
using Accessors
using Combinatorics
using SymmetryReduceBZ

include("crystal.jl")
include("interface.jl")
include("tb.jl")
include("nfe.jl")
include("misc.jl")

# Unit definitons

a₀ = 1.889726125Å
Ha = 27.2eV

export a₀, Ha, Å, eV
export Lattice, UnitCell, Crystal
export shiftenergy, fermilevel
export plotcrystal!
export evals, evecs, kvecs
export solve, plotSolution
export ReciprocalPath, DOS
export Hoppings, addhop!, addonsite!, addoverlap!
export TightBindingHamiltonian, PseudoPotentialHamiltonian
end
