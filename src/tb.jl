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

mutable struct Hoppings
    c::Crystal
    γs::Array{Hop}
    μs::Array{Onsite}
    Ss::Array{Overlap}
    maxij::Int

"""
    Hoppings(c::Crystal)

Create an empty hopping list based on crystal structure `c`.
"""
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

"""
    tbH(k, hops::Hoppings)

Returns a tuple `(H, S)` where `H` is the tight binding Hamiltonian matrix and `S` the overlap matrix at given momentum `k`.
"""
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
struct TightBindingProblem
    hops::Hoppings
    ks::Vector{Union{<:Real, Vector{<:Real}}}
    kp::Vector{<:Real}
    kl::Vector{Pair{Symbol, <:Real}}
end

"""
    TightBindingProblem(hops::Hoppings, kpositions::Vector, kstep::Real)

Assembles a tight binding problem that can be solved by the `solve` routine.

# Arguments
- `hops::Hoppings`: A hopping list.
- `kpositions::Vector{Pair(Symbol, T)}`: List of k-space branch vertices.
- `kstep::Real`: Step size along the k-space trajectory.
"""
function TightBindingProblem(hops::Hoppings, kpositions::Vector, kstep::Real)
    k = kpath(kpositions, kstep)
    TightBindingProblem(hops, k.path, k.plength, k.ppoints)
end

function solve(p::TightBindingProblem)::Solution
    h(k) = tbH(k, p.hops)
    sols = h.(p.ks) .|> x->eigen(x...)
    Solution(
        p.ks,
        p.kp,
        p.kl,
        [eig.values for eig ∈ sols],
        [mapslices(v->v./norm(v), eig.vectors, dims=1) for eig ∈ sols]
    )
end
