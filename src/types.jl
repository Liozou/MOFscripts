## Type definitions for intermediate representations of crystals up to a net

include("specialsolver.jl")

import Base: ==
using Serialization
using Tokenize

"""
    EquivalentPosition

Representation of a symmetry operation in 3D, defined by an affine function.
"""
struct EquivalentPosition
    mat::SMatrix{3,3,Rational{Int}, 9}
    ofs::SVector{3,Rational{Int}}
end

"""
Find the reference identifiers for the three dimensions for the CIF group called
symmetry_equiv_pos_as_xyz or space_group_symop_operation_xyz.
Usually this is simply ("x", "y", "z").
"""
function find_refid(eqs)
    isempty(eqs) && return ("x", "y", "z")
    for eq in eqs # Normally it should be the first one but I'm not sure of the specs here
        eq = filter(x -> !isspace(x), eq) # remove whitespaces
        ts = collect(tokenize(eq))
        any(Tokenize.Tokens.isoperator∘Tokenize.Tokens.kind, ts) && continue
        refid = [""]
        not_at_the_end = true
        for t in ts
            @assert not_at_the_end
            k = Tokenize.Tokens.kind(t)
            if k === Tokenize.Tokens.IDENTIFIER
                @assert refid[end] == ""
                refid[end] = Tokenize.Tokens.untokenize(t)
            elseif k === Tokenize.Tokens.COMMA || k === Tokenize.Tokens.SEMICOLON
                @assert refid[end] != ""
                push!(refid, "")
            elseif k === Tokenize.Tokens.ENDMARKER
                not_at_the_end = false
            else
                error("Input string {$eq} is not a valid symmetry equivalent")
            end
        end
        if not_at_the_end
            error("Unknown end of line marker for symmetry equivalent {$eq}")
        end
        if length(refid) != 3 || refid[end] == ""
            error("Input string {$eq} is not a valid symmetry equivalent")
        end
        return tuple(refid[1], refid[2], refid[3])
    end
    return ("x", "y", "z")
end

function Base.parse(::Type{EquivalentPosition}, s::AbstractString, refid=("x", "y", "z"))
    const_dict = Dict{String, Int}(refid[1]=>1, refid[2]=>2, refid[3]=>3)
    mat = zeros(Rational{Int}, 3, 3)
    ofs = zeros(Rational{Int}, 3)
    curr_num::Union{Int, Nothing} = nothing
    curr_val::Union{Rational{Int}, Nothing} = nothing
    curr_sign::Union{Bool, Nothing} = nothing
    encountered_div::Bool = false
    i = 1
    something_written = false
    for x in tokenize(s)
        k = Tokenize.Tokens.kind(x)
        k === Tokenize.Tokens.WHITESPACE && continue
        if k === Tokenize.Tokens.INTEGER
            @assert isnothing(curr_val)
            if encountered_div
                curr_val = Rational{Int}(Int(curr_num), parse(Int, x.val))
                curr_num = nothing
                encountered_div = false
            else
                @assert isnothing(curr_num)
                curr_num = parse(Int, x.val)
            end
        else
            @assert !encountered_div
            if k === Tokenize.Tokens.IDENTIFIER
                if !isnothing(curr_num)
                    @assert isnothing(curr_val)
                    curr_val = curr_num
                    curr_num = nothing
                end
                sign = isnothing(curr_sign) ? 1 : 2*curr_sign - 1
                val = isnothing(curr_val)  ? Rational{Int}(1) : curr_val
                j = const_dict[Tokenize.Tokens.untokenize(x)]
                mat[i,j] += sign * val
                curr_val = nothing
                curr_sign = nothing
                something_written = true
            else
                if x.kind === Tokenize.Tokens.FWD_SLASH
                    @assert isnothing(curr_val)
                    @assert !isnothing(curr_num)
                    encountered_div = true
                    continue
                end
                if !isnothing(curr_num)
                    @assert isnothing(curr_val)
                    curr_val = curr_num
                    curr_num = nothing
                end
                if !isnothing(curr_val)
                    sign = isnothing(curr_sign) ? 1 : 2*curr_sign - 1
                    @ifwarn begin if !iszero(ofs[i])
                            @warn "Existing offset already existing for position $i in {$s}"
                        end
                    end
                    ofs[i] += sign * Rational{Int}(curr_val)
                    curr_val = nothing
                    curr_sign = nothing
                else
                    @assert isnothing(curr_sign)
                end
                if x.kind === Tokenize.Tokens.PLUS
                    curr_sign = true
                elseif x.kind === Tokenize.Tokens.MINUS
                    curr_sign = false
                elseif k === Tokenize.Tokens.COMMA || k === Tokenize.Tokens.SEMICOLON
                    i > 2 && error("Too many dimensions specified for symmetry equivalent {$s}")
                    something_written || error("{$s} is not a symmetry equivalent (no dependency expressed in position $i)")
                    something_written = false
                    i += 1
                else
                    k !== Tokenize.Tokens.ENDMARKER && error("Unknown end of line marker for symmetry equivalent {$s}")
                    i!= 3 && error("Input string {$s} is not a valid symmetry equivalent")
                end
            end
        end
    end
    EquivalentPosition(SMatrix{3, 3, Rational{Int}, 9}(mat), SVector{3, Rational{Int}}(ofs))
end

function Base.show(io::IO, eq::EquivalentPosition)
    function rationaltostring(x::Rational{<:Integer}, notofs::Bool, first::Bool)
        if notofs && (x == 1 || x == -1)
            return x < 0 ? '-' : first ? "" : "+"
        end
        sign = x < 0 || first ? "" : "+"
        sign * (x.den == 1 ? string(x.num) : string(x.num)*'/'*string(x.den))
    end
    xyz = ('x', 'y', 'z')
    for i in 1:3
        first = true
        for j in 1:3
            if eq.mat[i,j] != 0
                coeff = rationaltostring(eq.mat[i,j], true, first)
                first = false
                print(io, coeff)
                print(io, xyz[j])
            end
        end
        if eq.ofs[i] != 0
            print(io, rationaltostring(eq.ofs[i], false, first))
        end
        i < 3 && print(io, ',')
    end
    nothing
end


"""
    Cell

Representation of a periodic cell in 3D. Contains information about the cell
(axes lengths and angles) and its symmetry group.
"""
struct Cell
    latticesystem::Symbol
    spacegroup::String
    tablenumber::Int
    mat::SMatrix{3,3,BigFloat,9} # Cartesian coordinates of a, b and c
    equivalents::Vector{EquivalentPosition}

    function Cell(lattice, space, table, a, b, c, α, β, γ, eq)
        cosα = cosd(α); cosβ = cosd(β); cosγ = cosd(γ); sinγ = sind(γ)
        ω = sqrt(1 - cosα^2 - cosβ^2 - cosγ^2 + 2*cosα*cosβ*cosγ)
        mat = SMatrix{3,3,BigFloat,9}([a  b*cosγ  c*cosβ ;
                                      0   b*sinγ  c*(cosα - cosβ*cosγ)/sinγ ;
                                      0   0       c*ω/sinγ ])
        return new(lattice, space, table, mat, eq)
    end

    function Cell(c::Cell, eqs::Vector{EquivalentPosition})
        return new(c.latticesystem, c.spacegroup, c.tablenumber, c.mat, eqs)
    end
end

function ==(c1::Cell, c2::Cell)
    c1.latticesystem == c2.latticesystem && c1.spacegroup == c2.spacegroup &&
    c1.tablenumber == c2.tablenumber && c1.mat == c2.mat && c1.equivalents == c2.equivalents
end

Cell() = Cell(Symbol(""), "P 1", 0, 10, 10, 10, 90, 90, 90, EquivalentPosition[])
function cell_parameters(mat::StaticArray{Tuple{3,3},BigFloat})
    _a, _b, _c = eachcol(mat)
    a = norm(_a)
    b = norm(_b)
    c = norm(_c)
    α = acosd(_b'_c/(b*c))
    β = acosd(_c'_a/(c*a))
    γ = acosd(_a'_b/(a*b))
    return (a, b, c, α, β, γ)
end
cell_parameters(cell::Cell) = cell_parameters(cell.mat)
function Cell(cell::Cell, mat::StaticArray{Tuple{3,3},BigFloat})
    a, b, c, α, β, γ = cell_parameters(mat)
    return Cell(cell.latticesystem, cell.spacegroup, cell.tablenumber, a, b, c, α, β, γ, cell.equivalents)
end
function Base.show(io::IO, cell::Cell)
    a, b, c, α, β, γ = Float64.(cell_parameters(cell))
    print(io, Cell, "(\"$(cell.spacegroup)\", ($a, $b, $c), ($α, $β, $γ))")
end


"""
    CIF

Representation of a .cif file.
"""
struct CIF
    cifinfo::Dict{String, Union{String, Vector{String}}}
    cell::Cell
    ids::Vector{Int}
    types::Vector{Symbol}
    pos::Matrix{Float64}
    bonds::Matrix{Bool}
end


function keep_atoms(cif::CIF, kept)
    kept_ids = sort!([cif.ids[i] for i in kept])
    unique!(kept_ids)
    idmap = Vector{Int}(undef, length(cif.types)) # upper bound on maximum(kept_ids)
    for (i,x) in enumerate(kept_ids)
        idmap[x] = i
    end
    return CIF(cif.cifinfo, cif.cell, [idmap[cif.ids[i]] for i in kept],
               cif.types[kept_ids], cif.pos[:, kept], cif.bonds[kept, kept])
end


function periodic_distance(u, v)
    dst = 0.0
    #=@inbounds=# for i in 1:3
        x = abs2(u[i] - v[i])
        if x > 0.25
            x = (1 - sqrt(x))^2
        end
        dst += x
    end
    return sqrt(dst)
end

"""
    periodic_distance(u, v, mat)

Distance between points u and v, given as triplets of fractional coordinates, in
a repeating unit cell of matrix mat.
The distance is the shortest between all equivalents of u and v.
"""
function periodic_distance(u, v, mat)
    dst = 0.0
    x = similar(u)
    uu = copy(u)
    vv = copy(v)
    #=@inbounds=# for i in 1:3
        diff = u[i] - v[i]
        if diff > 0.5
            x[i] = diff - 1
        elseif - diff > 0.5
            x[i] = diff + 1
        else
            x[i] = diff
        end
    end
    return norm(mat*x)
end


"""
    remove_partial_occupancy(::CIF)

Only keep one atom per atom site.
"""
function remove_partial_occupancy(cif::CIF)
    points::Vector{SVector{3,Float64}} = collect(eachcol(cif.pos))
    perm = sortperm(points)
    n = length(perm)
    last_triplet = points[perm[1]] .+ (1,1,1)
    last_position = 0
    toremove = Int[]
    for i in 1:n
        j = perm[i]
        thispoint = points[j]
        if norm(thispoint .- last_triplet) < 4e-4
            push!(toremove, max(last_position, j)) # keep the first that appears
        else
            last_position = j
            last_triplet = thispoint
        end
    end
    if isempty(toremove) # Nothing to do, we simply alias the CIF to help the compiler (it may be useless)
        return CIF(cif.cifinfo, cif.cell, cif.ids, cif.types, cif.pos, cif.bonds)
    end
    @ifwarn @warn "This CIF file contains a site with multiple atoms. Only one atom will be kept per atom site."
    sort!(toremove)
    ids = deleteat!(copy(cif.ids), toremove)
    tokeep = deleteat!(collect(1:n), toremove)
    bonds = cif.bonds[tokeep, tokeep]
    pos = cif.pos[:,tokeep]
    @assert allunique(collect(eachcol(pos)))
    return CIF(cif.cifinfo, cif.cell, ids, cif.types, pos, bonds)
end

"""
    prune_collisions(::CIF)

Remove all atoms that are suspiciously close to one another. This arises when all
the possible positions of at atom are superposed in the CIF file, typically for
a solvent which should be disregarded anyway.
"""
function prune_collisions(cif::CIF)
    toremove = Int[]
    points::Vector{SVector{3,Float64}} = collect(eachcol(cif.pos))
    n = length(points)
    for i in 1:n, j in (i+1):n
        if periodic_distance(points[i], points[j], Float64.(cif.cell.mat)) < 0.55
            push!(toremove, i)
            push!(toremove, j)
        end
    end
    if isempty(toremove) # Nothing to do, we simply alias the CIF to help the compiler (it may be useless)
        return CIF(cif.cifinfo, cif.cell, cif.ids, cif.types, cif.pos, cif.bonds)
    end
    @ifwarn @warn "This CIF file contains multiple colliding atoms. All colliding atoms will be removed."
    unique!(sort!(toremove))
    ids = deleteat!(copy(cif.ids), toremove)
    tokeep = deleteat!(collect(1:n), toremove)
    bonds = cif.bonds[tokeep, tokeep]
    pos = cif.pos[:,tokeep]
    return CIF(cif.cifinfo, cif.cell, ids, cif.types, pos, bonds)
end

"""
    expand_symmetry(::CIF)

Applies all the symmetry operations listed in the CIF file to the atoms and the bonds.
"""
function expand_symmetry(c::CIF)
    cif::CIF = remove_partial_occupancy(c)
    n = length(cif.ids)
    oldbonds = [(i,j) for i in 1:n for j in i:n if cif.bonds[i,j]]
    newbonds = Set{Tuple{Int,Int}}(oldbonds)
    newids::Vector{Int} = copy(cif.ids)
    newpos::Vector{Vector{Float64}} = collect(eachcol(cif.pos))
    ret = Vector{Vector{Int}}
    #=@inbounds=# for equiv in cif.cell.equivalents
        image = zeros(Int, n)
        for i in 1:length(cif.ids)
            v = newpos[i]
            p = Vector(equiv.mat*v .+ equiv.ofs)
            @. p = p - floor(p)
            for j in 1:length(newpos)
                if periodic_distance(newpos[j], p, Float64.(cif.cell.mat)) < 0.5
                    image[i] = j
                    break
                end
            end
            if image[i] == 0
                push!(newpos, p)
                push!(newids, cif.ids[i])
                image[i] = length(newpos)
            end
        end
        for (i,j) in oldbonds
            push!(newbonds, minmax(image[i], image[j]))
        end
    end
    m = length(newids)
    bonds = zeros(Bool, m, m)
    for (i,j) in newbonds
        bonds[i,j] = true
        bonds[j,i] = true
    end

    newcif = CIF(cif.cifinfo, deepcopy(cif.cell), newids, copy(cif.types), reduce(hcat, newpos), bonds)

    return prune_collisions(newcif)
end

# function set_unique_bond_type!(cif::CIF, bond_length, bonded_atoms::Tuple{Symbol, Symbol}, onlykeep, tol)
#     @assert iszero(cif.bonds)
#     indices = [i for i in 1:length(cif.ids) if cif.types[cif.ids[i]] ∈ onlykeep]
#     #=@inbounds=# for _i in 1:length(indices)
#     i = indices[_i]
#     cif.bonds[i,i] = false
#         Threads.@threads for _j in i+1:length(indices)
#             j = indices[_j]
#             if minmax(cif.types[cif.ids[i]], cif.types[cif.ids[j]]) == bonded_atoms
#                 bonded = abs2(periodic_distance(cif.pos[:,i], cif.pos[:,j], cif.cell.mat) - bond_length) <= tol
#                 cif.bonds[i,j] = bonded
#                 cif.bonds[j,i] = bonded
#             end
#         end
#     end
#     nothing
# end

"""
    edges_from_bonds(bonds, mat, pos)

Given an n×n adjacency matrix `bonds`, the 3×3 matrix of the cell `mat` and the
3×n matrix `pos` whose columns correspond to the positions of the atoms, extract
the list of PeriodicEdge3D corresponding to the bonds.
Since the adjacency matrix wraps bonds across the boundaries of the cell, the edges
are extracted so that the closest representatives are chosen to form bonds.
"""
function edges_from_bonds(bonds, mat, pos)
    @assert size(pos)[1] == 3
    n = size(pos)[2]
    edges = PeriodicEdge3D[]
    ref_dst = norm(mat*[1, 1, 1])
    for i in 1:n, k in (i+1):n
        bonds[k,i] || continue
        offset::Vector{SVector{3, Int}} = []
        old_dst = ref_dst
        for ofsx in -1:1, ofsy in -1:1, ofsz in -1:1
            # dst = norm(mat*(pos[:,i] .- (pos[:,k] .+ [ofsx, ofsy, ofsz])))
            dst = norm(pos[:,i] .- (pos[:,k] .+ (mat * [ofsx, ofsy, ofsz])))
            if abs2(dst - old_dst) < 1e-3
                push!(offset, (ofsx, ofsy, ofsz))
                old_dst = (dst + old_dst)/2
            elseif dst < old_dst
                empty!(offset)
                push!(offset, (ofsx, ofsy, ofsz))
                old_dst = dst
            end
        end
        for ofs in offset
            push!(edges, (i, k, ofs))
        end
    end
    return edges
end


"""
    Clusters

Classification of the atoms of a crystalline framework in different clusters.
For simple crystals, every atom is its own cluster.
For a MOF, a cluster is a SBU, which can be either organic or inorganic.
"""
struct Clusters
    sbus::Vector{Vector{PeriodicVertex3D}}
    classes::Vector{Int}
    attributions::Vector{Int}
    offsets::Vector{SVector{3,Int}}
end

function Clusters(n)
    sbus = [[PeriodicVertex3D(i)] for i in 1:n]
    classes = collect(1:n)
    attributions = collect(1:n)
    offsets = [zero(SVector{3,Int}) for _ in 1:n]
    return Clusters(sbus, classes, attributions, offsets)
end

Base.isempty(c::Clusters) = c.attributions == 1:length(c.attributions)

"""
    Crystal

Intermediate representation of a crystal, retaining information on the cell, and the
exact placement of the atoms and their type, as well as the residues which will be used as
vertices for the computation of the underlying topology.
"""
struct Crystal{T<:Union{Nothing,Clusters}}
    cell::Cell
    types::Vector{Symbol}
    clusters::T
    pos::Matrix{Float64}
    graph::PeriodicGraph3D
end

function ==(c1::Crystal{T}, c2::Crystal{T}) where T
    c1.cell == c2.cell && c1.types == c2.types && c1.clusters == c2.clusters &&
    c1.pos == c2.pos && c1.graph == c2.graph
end

Crystal{Nothing}(c::Crystal{Nothing}) = c
function Crystal{Nothing}(c::Crystal{Clusters})
    Crystal{Nothing}(c.cell, c.types, nothing, c.pos, c.graph)
end
function Crystal{Clusters}(c::Crystal, clusters::Clusters)
    Crystal{Clusters}(c.cell, c.types, clusters, c.pos, c.graph)
end

# function Crystal(cif::CIF)
#     c::CIF = expand_symmetry(cif)
#     if iszero(cif.bonds)
#         # The CIF files does not contain any bond, so we have to guess
#         atoms::Vector{Symbol} = unique!(sort(c1[].types))
#         if atoms == [:O, :Si] || atoms == [:Si] # zeolite
#             @info "Interpreting .cif file as a zeolite"
#             set_unique_bond_type!(c, 3.1, (:Si, :Si), (:Si,), 0.1)
#         elseif atoms == [:C] || atoms == [:C, :H]
#             set_unique_bond_type!(c, 1.54, (:C, :C), (:C,), 0.3)
#         elseif atoms == [:Cl, :Na] # toy example
#             set_unique_bond_type!(c, 2.82, (:Cl, :Na), (:Cl, :Na), 0.2)
#         else
#             error("Missing bonds on CIF object $(cif.cifinfo["data"])")
#         end
#     end
#     @assert !iszero(c.bonds)
#     bondedatoms::Vector{Int} = [i for i in 1:length(c.ids) if count(c.bonds[:,i]) > 1]
#     c = keep_atoms(c, bondedatoms) # Only keep atoms that have at least one bond
#     n = length(c.ids)
#     clusters = Clusters(n)
#     edges = edges_from_bonds(c.bonds, c.cell.mat, c.pos)
#     @assert !isempty(edges)
#     Crystal(c.cell, [c.types[x] for x in c.ids], clusters, c.pos, PeriodicGraph3D(n, edges))
# end


"""
    equilibrium(g::PeriodicGraph)

Return an equilibrium placement for the vertices of the graph, defined as a set
of positions such that each vertex is at the barycentre of its neighbors.
The returned equilibrium placement is such that the first vertex of the graph
is at the origin of the space.
"""
function equilibrium(g::PeriodicGraph{N}) where N
    n = nv(g)
    iszero(n) && return Matrix{Rational{Int128}}(undef, N, 0)
    Y = Matrix{Int}(undef, n, N)
    A = spzeros(Int, n, n)
    neigh = Vector{Int}(undef, n)
    offset = SizedVector{N,Int}(undef)
    for i in 1:n
        neigh .= 0
        offset .= 0
        count = 0
        for k in neighbors(g, i)
            k.v == i && continue
            count += 1
            neigh[k.v] += 1
            offset .-= k.ofs
        end
        Y[i,:] .= offset
        A[i,:] .= neigh
        A[i,i] = -count
    end

    return dixon_solve(Val(N), A[2:end,2:end], Y[2:end,:])
end

"""
    trim_topology(graph::PeriodicGraph)

Return a pair `(vmap, newgraph)` extracted from the input by removing vertices
of valence lower or equal to 1, and by replacing vertices of valence 2 by edges.

`vmap` maps the vertices of `newgraph` to their counterpart in `graph`.
"""
function trim_topology(graph::PeriodicGraph{N}) where N
    newgraph = PeriodicGraph{N}(nv(graph), collect(edges(graph)))
    remove_idx = Int[]
    flag = any(<=(2), degree(newgraph))
    vmap = collect(1:nv(graph))
    while flag # we alternate cycles when we remove valence 1 and valence 2
        # until no such vertex remains
        flag = any(<=(1), degree(newgraph))
        while flag
            n = nv(newgraph)
            for i in 1:n
                if degree(newgraph, i) <= 1
                    push!(remove_idx, i)
                end
            end
            map = rem_vertices!(newgraph, remove_idx)
            vmap = vmap[map]
            #pos = pos[vmap]
            empty!(remove_idx)
            flag = any(isone, degree(newgraph))
        end
        flag = any(==(2), degree(newgraph))

        while flag
            n = nv(newgraph)
            for i in 1:n
                if degree(newgraph, i) == 2
                    neigh1, neigh2 = neighbors(newgraph, i)
                    add_edge!(newgraph, PeriodicEdge{N}(neigh1.v, neigh2.v, neigh2.ofs .- neigh1.ofs))
                    push!(remove_idx, i)
                end
            end
            map = rem_vertices!(newgraph, remove_idx)
            vmap = vmap[map]
            #pos = pos[vmap]
            empty!(remove_idx)
            flag = any(==(2), degree(newgraph))
        end
        flag = any(<=(2), degree(newgraph))
    end
    return vmap, newgraph
end

trimmed_crystal(c::Crystal{Clusters}) = trimmed_crystal(coalesce_sbus(c))
function trimmed_crystal(c::Crystal{Nothing})
    vmap, graph = trim_topology(c.graph)
    return Crystal(c.cell, c.types[vmap], nothing, c.pos[:,vmap], graph)
end

"""
    CrystalNet{T<:Real}

Representation of a net as a topological abstraction of a crystal.

`T` is the numeric type used to store the exact coordinates of each vertex at the
equilibrium placement.
"""
struct CrystalNet{T<:Real}
    cell::Cell
    types::Vector{Symbol}
    pos::Vector{SVector{3,T}}
    graph::PeriodicGraph3D
end

function CrystalNet{T}(cell::Cell, types::AbstractVector{Symbol}, graph::PeriodicGraph3D,
                       placement::AbstractMatrix{T}) where T<:Real
    n = nv(graph)
    pos = Vector{SVector{3,T}}(undef, n)
    offsets = Vector{SVector{3,Int}}(undef, n)
    for (i, x) in enumerate(eachcol(placement))
        offsets[i] = floor.(Int, x)
        pos[i] = x .- offsets[i]
    end
    s = sortperm(pos)
    pos = pos[s]
    types = Symbol[types[s[i]] for i in 1:n]
    cell = Cell(cell, EquivalentPosition[])
    graph = offset_representatives!(graph, .-offsets)[s]
    # @assert all(pos[i] == mean(pos[x.v] .+ x.ofs for x in neighbors(graph, i)) for i in 1:length(pos))
    return CrystalNet(cell, types, pos, graph)
end

macro tryinttype(inttype, m, M, cell, types, graph, placement)
    quote
        if ((typemin($(esc(inttype))) <= $(esc(m))) & ($(esc(M)) <= typemax($(esc(inttype)))))
            return CrystalNet{Rational{$(esc(inttype))}}(
                        $(esc(cell)), $(esc(types)), $(esc(graph)),
                        Rational{$(esc(inttype))}.($(esc(placement))))
        end
    end
end

function trim_crystalnet!(graph, types, l, keep)
    tohandle::Vector{Int} = reduce(vcat, l)
    sort!(tohandle)
    toremove = keep ? deleteat!(collect(1:length(types)), tohandle) : tohandle
    vmap = rem_vertices!(graph, toremove)
    return vmap
    nothing
end

function CrystalNet(cell::Cell, types::AbstractVector{Symbol}, graph::PeriodicGraph3D)
    vmap, graph = trim_topology(graph)
    types = types[vmap]
    dimensions = PeriodicGraphs.dimensionality(graph)
    if haskey(dimensions, 0)
        @ifwarn @warn "Removing complex structure of dimension 0, possibly solvent residues."
        vmap = trim_crystalnet!(graph, types, dimensions[0], false)
        types = types[vmap]
        dimensions = PeriodicGraphs.dimensionality(graph)
    end

    if isempty(dimensions)
        return CrystalNet{Rational{Bool}}(cell, types, graph, Matrix{Rational{Bool}}(undef, 3, 0))
    end
    if !haskey(dimensions, 3)
        error("The input does not contain a 3D crystal.")
    end
    if length(dimensions) > 1 || only(keys(dimensions)) != 3
        @ifwarn @warn "Presence of periodic structures with a periodicity different than 3. At the moment, all substructures with periodicity other than 3 will be ignored."
        vmap = trim_crystalnet!(graph, types, get(dimensions, 3, Vector{Int}[]), true)
        types = types[vmap]
        dimensions = PeriodicGraphs.dimensionality(graph)
    end
    @assert !isempty(get(dimensions, 3, Vector{Int}[]))
    dim3 = dimensions[3]
    if length(dim3) > 1
        error("Multiple entertwined tridimensional structures. This is not currently handled.")
    end
    @assert sort!(only(dim3)) == collect(1:nv(graph))
    placement = equilibrium(graph)
    m = min(minimum(numerator.(placement)), minimum(denominator.(placement)))
    M = max(maximum(numerator.(placement)), maximum(denominator.(placement)))
    @tryinttype(Int8,   m, M, cell, types, graph, placement)
    @tryinttype(Int16,  m, M, cell, types, graph, placement)
    @tryinttype(Int32,  m, M, cell, types, graph, placement)
    @tryinttype(Int64,  m, M, cell, types, graph, placement)
    @tryinttype(Int128, m, M, cell, types, graph, placement)
    return CrystalNet{Rational{BigInt}}(cell, types, graph, placement)
    # Type-unstable function, but yields better performance than always falling back to Int128
end

function CrystalNet(c::Crystal{T}) where T
    if T === Clusters # the only alternative is c.clusters === nothing
        c = coalesce_sbus(c)
    end
    if ne(c.graph) == 0
        error("Empty graph. This probably means that the bonds have not been correctly attributed. Please switch to an input file containing explicit bonds.")
    end
    CrystalNet(c.cell, c.types, c.graph)
end

"""
    @enum ClusteringMode

Selection mode for the clustering of vertices. The choices are:
-   `InputClustering`: use the input residues as clusters. Fail if some atom does
    not belong to a residue.
-   `EachVertexClustering`: each vertex is its own cluster.
-   `MOFClustering`: discard the input residues and consider the input as a MOF. Identify
    organic and inorganic clusters using a simple heuristic based on the atom types.
-   `GuessClustering`: discard the input residues and try to identify the clusters as in `MOFClustering`.
    If it fails, use `EachVertexClustering`.
-   `AutomaticClustering`: if the input assigns each atom to a residue, use these residues
    as clusters. Otherwise, try to guess the clusters as in `GuessClustering`.
"""
@enum ClusteringMode begin
    InputClustering
    EachVertexClustering
    MOFClustering
    GuessClustering
    AutomaticClustering
end

function do_clustering(c::Crystal{T}, mode::ClusteringMode)::Tuple{Clusters, CrystalNet} where T
    n = length(c.types)
    if mode == InputClustering
        if T === Nothing
            throw(ArgumentError("Cannot use input residues as clusters: the input does not have residues."))
        else
            return Clusters(n), CrystalNet(c)
        end
    elseif mode == EachVertexClustering
        return Clusters(n), CrystalNet(Crystal{Nothing}(c))
    elseif mode == MOFClustering
        clusters = find_sbus(c)
        if length(clusters.sbus) <= 1
            throw(MissingAtomInformation("MOFClustering leads to a single cluster, choose a different clustering mode."))
        end
        return clusters, CrystalNet(Crystal{Clusters}(c, clusters))
    elseif mode == GuessClustering
        crystal = Crystal{Nothing}(c)
        try
            clusters, net = do_clustering(crystal, MOFClustering)
            if nv(net.graph) > 1
                return clusters, net
            end
        catch e
            if !(e isa MissingAtomInformation)
                rethrow()
            end
        end
        return do_clustering(crystal, EachVertexClustering)
    elseif mode == AutomaticClustering
        if T === Clusters
            return Clusters(n), CrystalNet(c)
        else
            return do_clustering(c, GuessClustering)
        end
    end
end

function CrystalNet(g::Union{PeriodicGraph3D,AbstractString,AbstractVector{PeriodicEdge3D}})
    graph = PeriodicGraph3D(g)
    cell = Cell()
    n = nv(g)
    types = [Symbol("") for _ in 1:n]
    return CrystalNet(cell, types, graph)
end


# function CrystalNet(cif::CIF)
#     @assert isempty(cif.cell.equivalents) # FIXME we only handle P1 representations
#     n = length(cif.ids)
#     @assert size(cif.pos)[1] == 3
#     @assert size(cif.pos)[2] == n
#
#     edges = edges_from_bonds(cif.bonds, cif.cell.mat, cif.pos)
#     pos = SVector{3,Rational{Int}}[rationalize.(Int,x) for x in eachcol(cif.pos)]
#     return CrystalNet(cif.cell, [cif.types[x] for x in cif.ids], pos, PeriodicGraph3D(edges))
# end
