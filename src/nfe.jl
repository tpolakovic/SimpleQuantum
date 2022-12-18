function nfH(k, n::Integer, V::Function, crystal::Crystal)
    G = crystal.lattice.G
    e = sort(Iterators.product(fill(-n:n,length(k))...) |> collect |> vec, by=norm) |> collect
    k = G' * k
    Gs = (G' * collect(g) for g ∈ e)
    H = [V(j - i) for i ∈ Gs, j ∈ Gs]
    H .+= diagm((1/2*norm(k .+ g)^2 for g ∈ Gs) |> collect)
end

function nfH(k, n::Integer, V::Matrix, crystal::Crystal)
    G = crystal.lattice.G
    e = sort(Iterators.product(fill(-n:n,length(k))...) |> collect |> vec, by=norm) |> collect
    k = G' * k
    Gs = (G' * g for g ∈ e)
    V .+ diagm((1/2*norm(k .+ g)^2 for g ∈ Gs) |> collect)
end

struct PseudoPotentialHamiltonian <: ReciprocalHamiltonian
    pot
    n::Integer
    c::Crystal
end

# interface

"""
     PsuedoPotentialHamiltonian(n::Integer, V::Function, c::Crystal)

Assembles the pseudopotential Hamiltonian.

The Hamiltonian is determined from the potential as a function of momentum `V` := V(k) with reciprocal lattice vectors from up to `n`-th shell Brillouin zone.

Can be callend on a vector in k-space to output the Hamiltonian.
"""
function PseudoPotentialHamiltonian(n::Integer, V::F, c::Crystal) where F <: Function
    PseudoPotentialHamiltonian(V, n, c)
end

"""
     PsuedoPotentialHamiltonian(V::Matrix, c::Crystal)

Assembles the pseudopotential Hamiltonian.

The potential terms of the Hamiltoniain are stored in the matrix `V`.

Can be callend on a vector in k-space to output the Hamiltonian.
"""
function PseudoPotentialHamiltonian(V::Matrix, c::Crystal)
    k = kpath(kpositions, kstep)
    (n1, n2) = size(V)
    if n1 != n2
        throw(ArgumentError("Pontential matrix is not square."))
    end
    n = n1
    PseudoPotentialHamiltonian(V, n, c)
end

# function (h::PseudoPotentialHamiltonian)(k::Vector)
#     nfH(k, h.n, h.pot, h.c)
# end

crystal(h::PseudoPotentialHamiltonian) = h.c

function getH(h::PseudoPotentialHamiltonian)
    k -> nfH(k, h.n, h.pot, h.c)
end

function Base.ndims(t::PseudoPotentialHamiltonian)
	t.c |> ndims
end
