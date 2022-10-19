"""
    Solution

Type containing data of the electronic structure calculation.
"""
struct Solution
    ks::Vector{Union{<:Real, Vector{<:Real}}}
    kp::Vector{<:Real}
    kl::Vector{Pair{Symbol, <:Real}}
    es::Vector{Vector{Real}}
    evs::Vector{Matrix{Complex}}
end

"""
    evals(s::Solution)

Returns an array of eigenvalues along the k-path.

The n-th element of the array contains an array of eigenvalues corresponding to n-th k-position returned by `kvecs`.
"""
evals(s::Solution) = s.es
"""
    evecs(s::Solution)

Returns an array of eigenvectors along the k-path.

Each element or the returned array is a matrix with eigenvectors as columns. n-th matrix corresponds to n-th k-position returned by `kvecs`.
"""
evecs(s::Solution) = s.evs
"""
    kvecs(s::Solution)

Returns the points along the k-space trajectory along which the solution was calculated.
"""
kvecs(s::Solution) = s.ks

"""
    solve(p)

Solves an assembled problem defined by constructors with name `<Model Name>Problem`.
"""
function solve(p)
    t = typeof(p)
    ArgumentError("$t has no solution implemented.")
end

function plotSolution(s::Solution)
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
