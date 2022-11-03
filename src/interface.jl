struct ReciprocalPath
    ks::Vector{Union{<:Real, Vector{<:Real}}}
    kp::Vector{<:Real}
    kl::Vector{Pair{Symbol, <:Real}}
end

abstract type ReciprocalHamiltonian end

function Base.ndims(t::T) where {T <: ReciprocalHamiltonian}
    e = typeof(t)
    throw(ArgumentError("$e needs to have `ndims` defined."))
end

function getH(t::T) where {T <: ReciprocalHamiltonian}
    e = typeof(t)
    throw(ArgumentError("$e needs to have `getH` defined."))
end

struct ReciprocalBandProblem
    h::ReciprocalHamiltonian
    kp::ReciprocalPath
end

function ReciprocalPath(kpositions::Vector, kstep::Real)
    k = kpath(kpositions, kstep)
    ReciprocalPath(k.path, k.plength, k.ppoints)
end

function (h::ReciprocalHamiltonian)(k)
    getH(h)(k)
end

function (rp::ReciprocalPath)(h::ReciprocalHamiltonian)
    ReciprocalBandProblem(h, rp)
end

function (h::ReciprocalHamiltonian)(rp::ReciprocalPath)
    ReciprocalBandProblem(h, rp)
end

# abstract type Solution end

"""
    BandSolution

Type containing data of the electronic band structure calculation.
"""
struct BandSolution #<: Solution
    ks::Vector{Union{<:Real, Vector{<:Real}}}
    kp::Vector{<:Real}
    kl::Vector{Pair{Symbol, <:Real}}
    evals::Vector{Vector{Real}}
    evecs::Vector{Matrix{Complex}}
end

"""
    evals(s::BandSolution)

Returns an array of eigenvalues of solution `s`.
"""
evals(s::BandSolution) = s.evals

"""
    evecs(s::BandSolution)

Returns an array of eigenvectors of solution `s`.
"""
evecs(s::BandSolution) = s.evecs

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
    for n ∈ 1:length(evals(s)[1])
        lines!(s.kp, [e[n] for e ∈ evals(s)])
    end
    fig[1,1] = ax
    fig
end

struct EFermi
    esincell
end

function(efs::EFermi)(s)
    es = evals(s)
    nks = es |> length
    esflat = es |> Iterators.flatten |> collect |> sort
    n = efs.esincell / 2
    n1 = floor(Int, n * nks)
    n2 = ceil(Int, n * nks)
    e = (esflat[n1] + esflat[n2]) / 2
    e * Ha
end

"""
    fermilevel(n)

Finds the Fermi level with `n` electrons in the unit cell.
"""
function fermilevel(n)
	EFermi(n)
end

"""
    shiftenergy(shift)

Rigidly shifts the energy axis by `shift`.
"""
function shiftenergy(shift)
	s -> @set s.evals = [e .- shift for e ∈ evals(s)]
end

function shiftenergy(shift::T) where T <: Quantity
	shift = Unitful.NoUnits(shift ./ Ha)
	shiftenergy(shift)
end

struct ReciprocalDOSProblem
	h
	dos
end

"""
    DOS(es, broadening, nq)

Sets up the density of states calculation on energies `es` with Lorentzian energy broadening `broadening` by solving the D-dimensional problem on an evenly spaced grid with `nq`^D points.
"""
struct DOS
    es
    broadening
    nq

    function DOS(es::Vector{<:Real}, broadening::Real, nq::Integer)
        new(es, broadening, nq)
    end
end

function DOS(es::Vector{<:Real}, nq::Integer)
    broadening = √(minimum(diff(es)))/50
    DOS(es, broadening, nq)
end

function DOS(es::Vector{<:Quantity}, broadening::T, nq::Integer) where T <: Quantity
    es = Unitful.NoUnits.(es ./ Ha)
    broadening = Unitful.NoUnits(broadening / Ha)
    DOS(es, broadening, nq)
end

function DOS(es::Vector{<:Quantity}, nq::Integer)
    es = Unitful.NoUnits.(es ./ Ha)
    DOS(es, nq)
end

function (d::DOS)(h::ReciprocalHamiltonian)
    ReciprocalDOSProblem(h, d)
end

struct DOSSolution
    es
    evals
    nks
    broadening
end

function (d::DOSSolution)(E)
    ϵ = d.broadening
    (1/(π * d.nks)) * map(d.evals) do ens
        map(ens) do En
            ϵ / ((E - En)^2 + ϵ^2)
        end |> sum
    end |> sum
end

evals(s::DOSSolution) = s.evals

function mpspace(q)
	[(2r - q - 1)/2q for r ∈ 1:q]
end

"""
    solve(dp::ReciprocalDOSProblem)

Calculates the density of states as defined in `dp`.
"""
function solve(dp::ReciprocalDOSProblem)
    h = dp.h
    dos = dp.dos
    q = dos.nq
    ϵ = dp.dos.broadening
    n = ndims(h)
    ks = reshape(Iterators.product(ntuple(_ -> mpspace(q), n)...) |> collect, 1, :)
    nks = length(ks)
    hs = h.(ks)

    if typeof(first(hs)) <: Tuple
        es = map(h -> eigen(h...).values, hs)
    else
        es = map(h -> eigen(h).values, hs)
    end

    DOSSolution(dos.es, es, nks, ϵ)
end

"""
    plotSolution(s::DOSSolution)

Plots the density of states of `s`.

The returned object is a figure that can be further modified if needed.
"""
function plotSolution(s::DOSSolution)
    fig = Figure()
    ax = Axis(fig)
    ax.xlabel = "E [Ha]"
    ax.ylabel = "DOS [a.u.]"

    lines!(s.es, s.(s.es))
    hidexdecorations!(ax, ticks=false, ticklabels=false, label=false)
    hideydecorations!(ax, label=false)

    fig[1,1] = ax
    fig
end
