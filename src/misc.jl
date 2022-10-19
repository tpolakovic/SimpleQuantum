"""
    kpath(kpoints, dk)

Constructs a k-space trajectory between multiple points.

Mostly used for finding branches between critical points in the Brillouin zone.

# Arguments
- `kpoints::Vector{Pair(:Symbol, Vector)}`: A pair of type <Position label> => <position in k-space>.
- `dk`: Spacing along the trajectory.
"""
function kpath(kpoints, dk)
    vertices = reverse([point.second for point ∈ kpoints])
    labels = [point.first for point ∈ kpoints]
    path = [last(vertices)]
    plength = zeros(typeof(dk),1)
    idxs = [1]

    while length(vertices) >= 2
        pop!(path)
        v1 = pop!(vertices)
        v2 = last(vertices)
        dir = v2 .- v1
        dirm = norm(dir)
        segment = [v1 .+ dir .* dd
            for dd ∈ range(start = 0, stop = 1, step = around(dk/dirm))]
        path = append!(path,segment)
        idxs = push!(idxs, last(idxs) + length(segment) - 1)
    end

    plength = append!(plength,
        Iterators.accumulate(
            +,[norm(v2 .- v1) for (v1,v2) ∈ zip(path[1:end-1], path[2:end])]
            ))
    points = [lab => plength[i] for (lab,i) ∈ zip(labels, idxs)]
    (path=path, plength=plength, ppoints=points)
end

function dedup_floats(itr)
    out = Vector{eltype(itr)}()
    push!(out, itr[1])
    for itrel ∈ itr
        if map(x -> !(itrel ≈ x), out) |> all
            push!(out, itrel)
        end
    end
    out
end

"""
    unique_neighbors(c::Crystal)

Finds unique nearest neighbor pairs in the crystal `c`.

This function *does not* return all nearest neighbor. If pair are conjugate (i.e. site1->site2 and site2->site1), only one will be returned.
"""
function unique_neighbors(c::Crystal)
    R = c.lattice.R
    positions = c.cell.positions
    supercellpositions = eltype(c.cell.positions)[]
    possible_hops = []
    out = []

    for (i, ipos) ∈ enumerate(positions)
        for offset ∈ Iterators.product(fill((-1,0,1), ndims(c.lattice))...)
            for (j, jpos) ∈ enumerate(positions)
                δ = (jpos .+ offset) - ipos
                    push!(possible_hops, (
                    i = i, j = j, δ = δ,
                    ))
            end
        end
    end
    filter!(v -> norm(v.δ) > 0, possible_hops)
    unique_dists = map(v -> (R * v.δ) |> norm, possible_hops) |> dedup_floats |> sort
    filter!(v -> norm(R * v.δ) <= unique_dists[1]+eps(), possible_hops)
    is_conjugate(u, v) = u.i == v.j && (u.δ .≈ -1 .* v.δ) |> all
    for hop ∈ possible_hops
        if map(v -> !is_conjugate(v, hop), out) |> all
            push!(out, hop)
        end
    end
    out
end
