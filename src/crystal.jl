# lattice
struct Lattice{T<:Real, N}
    R::Union{T, Array{T}}
    G::Union{T, Array{T}}
    V::T

    function Lattice(R::Union{<:Real,Matrix{<:Real}})
        G = 2π * inv(R)
        R,G = promote(R,G)
        V = det(R)
        dim = isempty(size(R)) ? 1 : first(size(R))
        new{eltype(R), dim}(R, G, V)
    end
end

function Lattice(R::Union{<:Quantity, Matrix{<:Quantity}})
    Unitful.NoUnits.(R ./ a₀) |> Lattice
end

"""
    Lattice(a, b, γ)

Construct a 2D Bravais lattice with lattice constants `a` and `b` and angle between them in degrees.

Lattice constants can be dimensionless or dimensionfull using `Unitful` units.
"""
function Lattice(a::T, b::T, γ::Real) where T <: Union{Real, Quantity}
    γ = deg2rad(γ)
    R = [a*sin(γ) zero(T);
         a*cos(γ) b]

    Lattice(R)
end

"""
    Lattice(a, b, c, α, β, γ)

Construct a 3D Bravais lattice with lattice constants `a`, `b` and `c` and angles `α`, `β` and `γ` between them in degrees.

Lattice constants can be dimensionless or dimensionfull using `Unitful` units.
"""
function Lattice(a::T, b::T, c::T, α::Real, β::Real, γ::Real) where T <: Union{Real, Quantity}
    α, β, γ = deg2rad.((α, β, γ))
    γ = (cos(α) * cos(β) - cos(γ)) / (sin(α) * sin(β))
    γ = clamp(γ, -1, 1) |> acos
    R = [a*sin(β) -b*sin(α)*cos(γ) zero(T);
         zero(T)  b*sin(α)*sin(γ)  zero(T);
         a*cos(β) b*cos(α)         c]

    Lattice(R)
end

Base.ndims(l::Lattice{T,N}) where {T,N} = N

# unit cell
struct UnitCell{T, N}
    positions::Vector{T}
    species::Vector{Symbol}

"""
    UnitCell(species::Vector{Symbol}, rs...)

Construct an crystal unit cell.

# Arguments
- `species::Vector{Symbol}`: Vector containing site identifiers as symbols. Must be the same length as number of site positions
- `rs::Union{Vector{T}, T}`: Site positions. Can be vectors or single numbers in 1D cells.

# Examples
```
UnitCell([:B, :N], [2//3, 1//3], [1//3, 2//3])
UnitCell{Vector{Rational{Int64}}, 2}(Vector{Rational{Int64}}[[2//3, 1//3], [1//3, 2//3]], [:B, :N])
```
"""
    function UnitCell(species::Vector{Symbol}, rs::T...) where T
        positions = collect(rs)
            if length(species) == length(rs)
                new{T, length(first(rs))}(positions, species)
            else
                throw(DimensionMismatch("Number of species and positions does not match"))
            end
    end

end

"""
    UnitCell(species::Vector{Symbol}, rs...)

Construct an crystal unit cell.

# Arguments
- `species::Symbol`: Site identifier as symbol. All sites will have this identifier.
- `rs::Union{Vector{T}, T}`: Site positions. Can be vectors or single numbers in 1D cells.

# Examples
```
UnitCell(:C, [2//3, 1//3], [1//3, 2//3])
UnitCell{Vector{Rational{Int64}}, 2}(Vector{Rational{Int64}}[[2//3, 1//3], [1//3, 2//3]], [:C, :C])
```
"""
function UnitCell(species::Symbol, rs...)
    species = fill(species, length(rs))
    UnitCell(species, rs...)
end

"""
    UnitCell(rs...)

Construct an crystal unit cell. Sites will have a generic label `:none`

# Arguments
- `rs::Union{Vector{T}, T}`: Site positions. Can be vectors or single numbers in 1D cells.

# Examples
```
UnitCell([2//3, 1//3], [1//3, 2//3])
UnitCell{Vector{Rational{Int64}}, 2}(Vector{Rational{Int64}}[[2//3, 1//3], [1//3, 2//3]], [:none, :none])
```
"""
function UnitCell(rs...)
    UnitCell(:none, rs...)
end

Base.ndims(c::UnitCell{T,N}) where {T,N} = N
Base.length(c::UnitCell) = c.positions |> length

# crystal
"""
    Crystal(L::Lattice, c::UnitCell)

Creates a crystal structure out of unit cell `c` on a Bravais lattice `L`.
"""
struct Crystal{N}
    lattice::Lattice
    cell::UnitCell

    function Crystal(l::Lattice{T,N}, c::UnitCell) where {T,N}
        new{N}(l, c)
    end
end

function plotcrystal!(ax, c::Crystal, vertpos; ncells, showcell, showbonds, cmap)
    R = c.lattice.R

    if showbonds
        for offset ∈ Iterators.product(fill(0:(ncells-1), ndims(c.lattice))...)
            supercellpositions = eltype(c.cell.positions)[]
            for pos ∈ c.cell.positions
                for tn ∈ Iterators.product(fill((-1,0,1), ndims(c.lattice))...)
                    push!(supercellpositions, R * (pos .+ tn))
                end
            end

            for pos ∈ c.cell.positions
                pos = R * pos
                rs = sort(map(v -> v .- pos, supercellpositions), by=norm)
                filter!(v -> norm(v) > 0, rs)
                nns = filter(v -> norm(v) ≈ norm(first(rs)), rs)

                for nn ∈ nns
                    pos1 = pos .+ R * collect(offset)
                    pos2 = pos .+ nn .+ R * collect(offset)
                    lines!(ax, [tuple(pos1...), tuple(pos2...)], color=:black, linewidth=2)
                end
            end
        end
    end

    for offset ∈ Iterators.product(fill(0:(ncells-1), ndims(c.lattice))...)
        sps = unique(c.cell.species)
        nsps = length(sps)
        cmap = Makie.categorical_colors(cmap, nsps > 1 ? nsps : 2)
        if showcell
            cellvertices = map(v -> tuple((R * (v .+ offset))...), eachcol(vertpos))
            lines!(ax, cellvertices, color=:grey, linewidth=0.5)
        end
        for (i,sp) ∈ enumerate(sps)
            idxs = findall(x -> x == sp, c.cell.species)
            pos = [tuple((R * (c.cell.positions[i] .+ offset))...) for i ∈ idxs]
            meshscatter!(ax, pos, markersize=0.2, color=cmap[i], shading=true, label=string(sp))
        end
    end

    nothing
end

Base.ndims(::Crystal{N}) where {N} = N

plotcrystal!(ax, c::Crystal{2}; ncells=1, showcell=true, showbonds=true, cmap=:PuOr_5) =
    plotcrystal!(ax, c, [0 0 1 1 0; 0 1 1 0 0]; ncells, showcell, showbonds, cmap=cmap)

plotcrystal!(ax, c::Crystal{3}; ncells=1, showcell=true, showbonds=true, cmap=:PuOr_5) = plotcrystal!(ax, c,
    [0 1 1 0 0 0 0 0 0 1 1 1 1 1 1 1 0 0 1;
     0 0 1 1 0 0 1 1 0 0 0 1 1 0 0 1 1 0 0;
     0 0 0 0 0 1 1 0 0 0 1 1 0 0 1 1 1 1 1]; ncells, showcell, showbonds, cmap=cmap)
