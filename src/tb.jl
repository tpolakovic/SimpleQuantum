struct Hop
    γ::Number
    i::Int
    j::Int
    offset::Union{<:Real, Array{<:Real}}
end

struct Onsite
    μ::Number
    ii::Int
end

struct Overlap
    S::Number
    i::Int
    j::Int
    offset::Union{<:Real, Array{<:Real}}
end

"""
    Hoppings(c::Crystal)

Create an empty hopping list based on crystal structure `c`.
"""
mutable struct Hoppings
    c::Crystal
    γs::Array{Hop}
    μs::Array{Onsite}
    Ss::Array{Overlap}
    maxij::Int

    function Hoppings(c::Crystal)
        new(c, [], [], [], 0)
    end
end

"""
    addhop!(hops::Hoppings, γ, i, j, δ)

Appends a hopping into `hops` with amplitude γ between sites `i` and `j` separated by `δ`.

The separation `δ` is in reduced (fractional) coordinates.

# Warning
- Hoppings are conjugated automatically. This means that you need to only provide hopping i->j and j->i will be generated for you.
"""
function addhop!(hops, γ, i, j, δ)
    γ = Unitful.NoUnits(γ / Ha)
    offset = δ
    push!(hops.γs, Hop(γ, i, j, offset))
    hops.maxij = max(i, j, hops.maxij)
    hops
end

"""
    addonsite!(hops::Hoppings, μ, ii)

Adds onsite energy `μ` to site `ii`.
"""
function addonsite!(hops, μ, ii)
    μ = Unitful.NoUnits(μ / Ha)
    push!(hops.μs, Onsite(μ, ii))
    hops.maxij = max(hops.maxij, ii)
    hops
end

"""
    addhop!(hops::Hoppings, γ, i, j, δ)

Appends an overlap into `hops` with value S between sites `i` and `j` separated by `δ`.

The separation `δ` is in reduced (fractional) coordinates.

"""
function addoverlap!(hops, S, i, j, δ)
    offset = δ
    push!(hops.Ss, Overlap(S, i, j, offset))
    hops.maxij = max(i, j, hops.maxij)
    hops
end

function tbH(k, hops::Hoppings)
    n = hops.maxij
    ham = zeros(Complex,n,n)
    Smat = zeros(ComplexF64,n,n)

    for hop ∈ hops.γs
        i, j = hop.i, hop.j
        r = hop.offset
        γ = hop.γ * exp(2im * π * (k ⋅ r))
        ham[i,j] += γ
        ham[j,i] += conj(γ)
    end

    for onsite ∈ hops.μs
        i = onsite.ii
        ham[i,i] = onsite.μ
    end

    for overlap ∈ hops.Ss
        i, j = overlap.i, overlap.j
        r = overlap.offset
        S = overlap.S * exp(2im * π * (k ⋅ r))
        Smat[i,j] += S
        Smat[j,i] += conj(S)
    end
    Smat += diagm(ones(n))
    (ham, Smat)
end

# interface
"""
    TightBindingProblem(hops::Hoppings)

Tight binding problem definition.

Can be called on a vector in k-space to output the Hamiltonian.
"""
struct TightBindingHamiltonian <: ReciprocalHamiltonian
    hops::Hoppings
end

# function (p::TightBindingHamiltonian)(k)
#     tbH(k, p.hops)
# end

function getH(h::TightBindingHamiltonian)
    k -> tbH(k, h.hops)
end

function Base.ndims(t::TightBindingHamiltonian)
	t.hops.c |> ndims
end
