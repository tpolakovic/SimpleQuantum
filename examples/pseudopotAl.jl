using SimpleQuantum
using LinearAlgebra
using GLMakie

# Define the aluminum crystal structure.
Al = Crystal(
    Lattice(2.856Å, 2.856Å, 2.856Å, 60, 60, 60),
    UnitCell(:Al, [0.0, 0.0, 0.0])
)

# Momentum representation of a Thomas-Fermi potential with charge of Q = 3 and screening length |q| = 3.
V(k) = ifelse(norm(k) ≈ 0, 0, 4π * 3/(norm(k)^2 .+ 3^2))

alH = PseudoPotentialHamiltonian(2, V, Al)

kpath = SimpleQuantum.ReciprocalPath([
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

alH |> kpath |> solve |> plotSolution

ylims!(current_axis(), (0,10))
