# SimpleQuantum

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
    addhop!(grhops, -1.0Ha, hop.i, hop.j, hop.δ)
end

# Assemble the tight binding problem.
grprob = TightBindingProblem(grhops, [
    :K => [1/3,1/3],
    :Γ => [0,0],
    :M => [1/2,0],
    :K => [1/3,1/3]
], 0.005)

# Solve the problem.
sol = solve(grprob)

# Plot the band diagram.
plotSolution(sol)
```
Output:

<img src="https://user-images.githubusercontent.com/9288586/196821411-d52fd1a1-5ca1-487a-8295-637948a6b750.png" width=50%>

### Nearly free electron model of aluminum

``` julia
using SimpleQuantum
using LinearAlgebra
using GLMakie

Al = Crystal(
    Lattice(2.856Å,2.856Å,2.856Å,60,60,60),
    UnitCell(:Al, [0.0,0.0,0.0])
)

# Momentum representation of a Thomas-Fermi potential with charge of Q = 3 and screening length |q| = 3.
V(k) = ifelse(norm(k) ≈ 0, 0, 4π * 3/(norm(k)^2 .+ 3^2))

# Use reciprocal vectors of magnitude up to 2 reduced momentum units.
alprob = NearlyFreeElectronProblem(2, V, Al, [
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

sol = solve(alprob)

plotSolution(sol)
ylims!(current_axis(), (0,10))
       
```

output:

<img src="https://user-images.githubusercontent.com/9288586/197043712-79040f5c-5b09-48c4-ab62-302871eb418f.png" width=50%>
