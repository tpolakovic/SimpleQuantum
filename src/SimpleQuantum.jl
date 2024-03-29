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
using MinkowskiReduction

include("crystal.jl")
include("interface.jl")
include("tb.jl")
include("nfe.jl")
include("misc.jl")

# Unit definitons

a₀ = 0.529177249Å
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
