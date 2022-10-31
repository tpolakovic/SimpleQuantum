struct ReciprocalPath
    ks::Vector{Union{<:Real, Vector{<:Real}}}
    kp::Vector{<:Real}
    kl::Vector{Pair{Symbol, <:Real}}
end

abstract type ReciprocalHamiltonian end

struct ReciprocalBandProblem
    h::ReciprocalHamiltonian
    kp::ReciprocalPath
end

function ReciprocalPath(kpositions::Vector, kstep::Real)
    k = kpath(kpositions, kstep)
    ReciprocalPath(k.path, k.plength, k.ppoints)
end

function (rp::ReciprocalPath)(h::ReciprocalHamiltonian)
    ReciprocalBandProblem(h, rp)
end

abstract type Solution end

"""
    evals(s::Solution)

Returns an array of eigenvalues of solution `s`.
"""
evals(s::Solution) = s.es

"""
    evecs(s::Solution)

Returns an array of eigenvectors of solution `s`.
"""
evecs(s::Solution) = s.evs

"""
    BandSolution

Type containing data of the electronic band structure calculation.
"""
struct BandSolution <: Solution
    ks::Vector{Union{<:Real, Vector{<:Real}}}
    kp::Vector{<:Real}
    kl::Vector{Pair{Symbol, <:Real}}
    es::Vector{Vector{Real}}
    evs::Vector{Matrix{Complex}}
end

"""
    kvecs(s::BandSolution)

Returns the points along the k-space trajectory along which the solution was calculated.
"""
kvecs(s::BandSolution) = s.ks

"""
    solve(rbp::ReciprocalBandProblem)

Calculates the electronic bands of `rbp`.
"""
function solve(rbp::ReciprocalBandProblem)
    kbranch = rbp.kp
    hs = rbp.h.(kbranch.ks)
    if typeof(first(hs)) <: Tuple
        sols = hs .|> x -> eigen(x...)
    else
        sols = hs .|> eigen
    end
    BandSolution(
        kbranch.ks,
        kbranch.kp,
        kbranch.kl,
        [eig.values for eig ∈ sols],
        [mapslices(v -> v ./ norm(v), eig.vectors, dims = 1) for eig ∈ sols]
    )
end

"""
    plotSolution(s::BandSolution)

Plots the band diagram of `s`.

The returned object is a figure that can be further modified if needed.
"""
function plotSolution(s::BandSolution)
    fig = Figure()
    ax = Axis(fig)
    ax.xticks = ([p.second for p ∈ s.kl],
                 [string(p.first) for p ∈ s.kl])
    ax.ylabel = "E [Ha]"

    xlims!(ax, (0, s.kp[end]))
    hideydecorations!(ax, ticks=false, ticklabels=false, label=false)
    for n ∈ 1:length(s.es[1])
        lines!(s.kp, [e[n] for e ∈ s.es])
    end
    fig[1,1] = ax
    fig
end
