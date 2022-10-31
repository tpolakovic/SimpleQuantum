using SimpleQuantum
using LinearAlgebra
using Match

# Define the orbitals at each site.
orbs = [
    (:s, [2/3, 1/3]),  #1
    (:s, [1/3, 2/3]),  #2
    ###
    (:px, [2/3, 1/3]), #3
    (:px, [1/3, 2/3]), #4
    ###
    (:py, [2/3, 1/3]), #5
    (:py, [1/3, 2/3]), #6
    ###
    (:pz, [2/3, 1/3]), #7
    (:pz, [1/3, 2/3]), #8
]

# Construct an all-electron graphene crystal.
grorbs = Crystal(
    Lattice(2.468Å, 2.468Å, 120),
    UnitCell([first(orb) for orb ∈ orbs], (last(orb) for orb ∈ orbs)...)
)

# Define the hopping/onsite energies and ovelap integrals.
hopping_ints = Dict(
    :ssσ => -5.729eV,
    :spσ => 5.618eV,
    :ppσ => 6.050eV,
    :ppπ => -3.070eV
)
onsite_es = Dict(
    :s => -8.37eV,
    :p => 0eV
)
overlap_vals = Dict(
    :ssσ => 0.102,
    :spσ => -0.171,
    :ppσ => -0.377,
    :ppπ => 0.070
)

# Initialize an empty hopping list based on the crystal structure.
grhops = Hoppings(grorbs)

# Fill the hopping list
dir_cos(r) = r ./ norm(r)
for hop ∈ SimpleQuantum.unique_neighbors(grorbs)
    orb_types = (orbs[hop.i][1], orbs[hop.j][1])
    l, m = dir_cos(grorbs.lattice.R * hop.δ)
    γ = @match orb_types begin
        (:s, :s)   => hopping_ints[:ssσ]
        (:s, :px)  => l * hopping_ints[:spσ]
        (:s, :py)  => m * hopping_ints[:spσ]
        #(:s, :pz) => 0.
        (:px, :s)  => l * hopping_ints[:spσ]
        (:px, :px) => l^2 * hopping_ints[:ppσ] + (1-l^2) * hopping_ints[:ppπ]
        (:px, :py) => l*m * (hopping_ints[:ppσ] - hopping_ints[:ppπ])
        #(:px, :pz) => 0.
        (:py, :s)  => m * hopping_ints[:spσ]
        (:py, :px) => l*m * (hopping_ints[:ppσ] - hopping_ints[:ppπ])
        (:py, :py) => m^2 * hopping_ints[:ppσ] + (1-m^2) * hopping_ints[:ppπ]
        #(:py, :pz) => 0.
        #(:pz, :s)  => 0.
        #(:pz, :px) => 0.
        #(:pz, :py) => 0.
        (:pz, :pz) => hopping_ints[:ppπ]
        _ => 0eV
    end
    addhop!(grhops, γ, hop.i, hop.j, hop.δ)

    s = @match orb_types begin
        (:s, :s)   => overlap_vals[:ssσ]
        (:s, :px)  => l * overlap_vals[:spσ]
        (:s, :py)  => m * overlap_vals[:spσ]
        #(:s, :pz) => 0.
        (:px, :s)  => l * overlap_vals[:spσ]
        (:px, :px) => l^2 * overlap_vals[:ppσ] + (1-l^2) * overlap_vals[:ppπ]
        (:px, :py) => l*m * (overlap_vals[:ppσ] - overlap_vals[:ppπ])
        #(:px, :pz) => 0.
        (:py, :s)  => m * overlap_vals[:spσ]
        (:py, :px) => l*m * (overlap_vals[:ppσ] - overlap_vals[:ppπ])
        (:py, :py) => m^2 * overlap_vals[:ppσ] + (1-m^2) * overlap_vals[:ppπ]
        #(:py, :pz) => 0.
        #(:pz, :s)  => 0.
        #(:pz, :px) => 0.
        #(:pz, :py) => 0.
        (:pz, :pz) => overlap_vals[:ppπ]
        _ => 0
    end
    addoverlap!(grhops, s, hop.i, hop.j, hop.δ)
end

for (i, orb) ∈ enumerate(orbs)
    orb_type = orb[1]
    μ = @match orb_type begin
        :s => onsite_es[:s]
        _  => onsite_es[:p]
    end
    addonsite!(grhops, μ , i)
end

# Define the tight binding Hamiltonian
grH = TightBindingHamiltonian(grhops)

# Define the momentum path.
kpath = SimpleQuantum.ReciprocalPath([
           :K => [1/3,1/3],
           :Γ => [0,0],
           :M => [1/2,0],
           :K => [1/3,1/3]
       ], 0.005)

# Solve and plot the problem.
grH |> kpath |> solve |> plotSolution
