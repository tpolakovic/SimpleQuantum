function nfH(k, n::Integer, V::Function, crystal::Crystal)
    G = crystal.lattice.G
    e = sort(Iterators.product(fill(-n:n,length(k))...) |> collect |> vec, by=norm) |> collect
    k = G' * k
    Gs = (G' * collect(g) for g ∈ e)
    H = [V(j - i) for i ∈ Gs, j ∈ Gs]
    H .+= diagm((1/2*norm(k .+ g)^2 for g ∈ Gs) |> collect)
end

function nfH(k, n::Integer, V::Matrix, crystal::Crystal)
    #(n, _) = size(V)
    G = crystal.lattice.G
    e = sort(Iterators.product(fill(-n:n,length(k))...) |> collect |> vec, by=norm) |> collect
    k = G' * k
    Gs = (G' * g for g ∈ e)
    V .+ diagm((1/2*norm(k .+ g)^2 for g ∈ Gs) |> collect)
end

struct NearlyFreeElectronProblem
    pot
    n::Integer
    c::Crystal
    ks::Vector{Union{<:Real, Vector{<:Real}}}
    kp::Vector{<:Real}
    kl::Vector{Pair{Symbol, <:Real}}
end

function NearlyFreeElectronProblem(n::Integer, V::F, c::Crystal, kpositions::Vector, kstep::Real) where F <: Function
    k = kpath(kpositions, kstep)
    NearlyFreeElectronProblem(V, n, c, k.path, k.plength, k.ppoints)
end

function NearlyFreeElectronProblem(V::Matrix, c::Crystal, kpositions::Vector, kstep::Real)
    k = kpath(kpositions, kstep)
    (n1, n2) = size(V)
    if n1 != n2
        throw(ArgumentError("Pontential matrix is not square."))
    end
    n = n1
    NearlyFreeElectronProblem(V, n, c, k.path, k.plength, k.ppoints)
end

function solve(p::NearlyFreeElectronProblem)::Solution
    h(k) = nfH(k, p.n, p.pot, p.c)
    sols = h.(p.ks) .|> eigen
    Solution(
        p.ks,
        p.kp,
        p.kl,
        [eig.values for eig ∈ sols],
        [mapslices(v->v./norm(v), eig.vectors, dims=1) for eig ∈ sols]
    )
end
