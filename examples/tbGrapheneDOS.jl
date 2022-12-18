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

# Spec up the DOS calculation: Energy in range (-0.5 Ha, 0.5 Ha) with step of 10 mHa and 151x151 grid points.
grdos = DOS([i * Ha for i ∈ -0.5:0.01:0.5], 151)

# Solve the problem and plot the DOS.
grH |> grdos |> solve |> plotSolution
