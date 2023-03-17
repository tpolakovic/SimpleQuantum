using SimpleQuantum
using LinearAlgebra
using GLMakie

# Define the aluminum crystal structure.
Al = Crystal(
    Lattice(2.856Å, 2.856Å, 2.856Å, 60, 60, 60),
    UnitCell(:Al, [0.0, 0.0, 0.0])
)

# Momentum representation of a Thomas-Fermi potential with charge of Q = 3 and screening length |q| = 3.
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

# Scale the energy axis to focus only on reasonable states.
ylims!(current_axis(), (-0.2,0.4))
