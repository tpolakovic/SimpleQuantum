# SimpleQuantum

A Julia package for calculating properties of crystalline solids. Currently supports tight binding calculations on up to 3-dimensional problems.

## Installation

In Julia REPL, press `]` to enter pkg mode and enter:

```
pkg> add https://github.com/tpolakovic/SimpleQuantum
```
## Sample

To calculate the two-band graphene structure:

``` julia

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
