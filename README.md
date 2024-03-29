# SimpleQuantum

[![DOI](https://zenodo.org/badge/554441370.svg)](https://zenodo.org/badge/latestdoi/554441370)

A Julia package for calculating properties of crystalline solids. Currently supports tight binding and nearly free electron model calculations on up to 3-dimensional problems.

Read the posts on the [Computational Physics for the Masses series](https://tpolakovic.github.io) for more detailed description of code and algorithms.

## Installation

In Julia REPL, press `]` to enter pkg mode and enter:

```
pkg> add https://github.com/tpolakovic/SimpleQuantum
```

## Sample

### Tight binding model of two-band graphene:

``` julia
using SimpleQuantum

# Define the graphene crystal.
graphene = Crystal(
    Lattice(2.468Å, 2.468Å, 120),
    UnitCell(:C, [2/3, 1/3], [1/3, 2/3])
)

# Create an empty hopping list based on graphene.
grhops = Hoppings(graphene)

# Iterate through unique nearest neighbor pairs and add them to the hopping list.
for hop ∈ SimpleQuantum.unique_neighbors(graphene)
    addhop!(grhops, -2.8eV, hop.i, hop.j, hop.δ)
end

# Define the tight binding Hamiltonian
grH = TightBindingHamiltonian(grhops)

# Define the momentum path.
kpath = ReciprocalPath([
           :K => [1/3,1/3],
           :Γ => [0,0],
           :M => [1/2,0],
           :K => [1/3,1/3]
       ], 0.005)

# Solve and plot the problem.
grH(kpath) |> solve |> plotSolution
```
Output:

<img src="https://user-images.githubusercontent.com/9288586/199577501-a192bdaf-3da7-4186-b749-2c56ef05cc9d.png" width=50%>

### Nearly free electron model of aluminum:

``` julia
using SimpleQuantum
using LinearAlgebra
using GLMakie

# Define the aluminum crystal structure.
Al = Crystal(
    Lattice(2.856Å, 2.856Å, 2.856Å, 60, 60, 60),
    UnitCell(:Al, [0.0, 0.0, 0.0])
)

# Momentum representation of a Thomas-Fermi potential with charge of Q = 3 and screening length |q| = 10.
V(k) = ifelse(norm(k) ≈ 0, 0, 4π * 3/(norm(k)^2 .+ 10^2))

# Define the Hamiltonian using the pseudopotential above with reciprocal lattice vectors from up to the 2nd shell
alH = PseudoPotentialHamiltonian(2, V, Al)

# Define the momentum path.
kpath = ReciprocalPath([
    :Γ => [0,0,0],
    :X => [1/2,0,1/2],
    :W => [1/2,1/4,3/4],
    :K => [3/8,3/8,3/4],
    :Γ => [0,0,0],
    :L => [1/2,1/2,1/2],
    :U => [5/8,1/4,5/8],
    :W => [1/2,1/4,3/4],
    :L => [1/2,1/2,1/2],
    :K => [3/8,3/8,3/4]
], 0.01)

# Solve the problem
sol = alH(kpath) |> solve

# Find the Fermi energy (3 free electrons in the unit cell of Al).
ef = sol |> fermilevel(3)
println("Ef = $ef")

# Plot the band diagram with energy shifted such that E_Fermi = 0.
sol |> shiftenergy(ef) |> plotSolution

ylims!(current_axis(), (-0.2,0.4))
```

output:

<img src="https://user-images.githubusercontent.com/9288586/225970346-c0e4b431-137c-4b7a-8d5a-459bbf216cdc.png" width=50%>
