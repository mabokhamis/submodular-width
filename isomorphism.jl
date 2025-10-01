module Isomorphism

using Combinatorics

using ..HypergraphWidths

export isomorphic_hash, are_isomorphic, isomorphically_unique

function _hash(x::Vector{Vector{UInt64}})::UInt64
    y = map(v -> Base.hash(sort(v)), x)
    return Base.hash(sort(y))
end

function isomorphic_hash(H::Hypergraph{T}) where T
    neighbors = Dict{T, Set{Set{T}}}()
    for v in H.vars
        neighbors_v = Set{Set{T}}()
        for e in H.edges
            if v in e
                push!(neighbors_v, Set{T}(u for u in e if u != v))
            end
        end
        neighbors[v] = neighbors_v
    end
    h = Dict{T, UInt64}(v => UInt64(0) for v in H.vars)
    num_distinct = 1
    while true
        new_h = Dict{T, UInt64}()
        for (v, neighbors_v) in neighbors
            h_neighbors_v = Vector{UInt64}[UInt64[h[u] for u in e] for e in neighbors_v]
            new_h[v] = Base.hash(h[v], _hash(h_neighbors_v))
        end
        new_num_distinct = length(Set{UInt64}(values(new_h)))
        new_num_distinct == num_distinct && break
        h = new_h
        num_distinct = new_num_distinct
    end
    return Base.hash(sort(collect(values(h)))), h
end

function are_isomorphic(H1::Hypergraph{T}, H2::Hypergraph{T}) where T
    h1, map1 = isomorphic_hash(H1)
    h2, map2 = isomorphic_hash(H2)
    h1 != h2 && return false
    rev1 = Dict{UInt64, Set{T}}()
    for (v, h) in map1
        @assert h isa UInt64
        if !haskey(rev1, h)
            rev1[h] = Set{T}()
        end
        push!(rev1[h], v)
    end
    rev2 = Dict{UInt64, Set{T}}()
    for (v, h) in map2
        @assert h isa UInt64
        if !haskey(rev2, h)
            rev2[h] = Set{T}()
        end
        push!(rev2[h], v)
    end
    Set{UInt64}(keys(rev1)) != Set{UInt64}(keys(rev2)) && return false
    vars1 = T[]
    vars2 = Vector{T}[]
    for k in keys(rev1)
        length(rev1[k]) != length(rev2[k]) && return false
        append!(vars1, collect(rev1[k]))
        push!(vars2, collect(rev2[k]))
    end
    perm2 = Iterators.product(map(vs -> permutations(vs), vars2)...)
    for p in perm2
        vars2 = Vector{T}(vcat(p...))
        _are_isomorphic(H1, H2, vars1, vars2) && return true
    end
    return false
end

function isomorphically_unique(H::Vector{Hypergraph{T}}) where T
    map = Dict{UInt64, Set{Hypergraph{T}}}()
    for G in H
        h = isomorphic_hash(G)[1]
        @assert h isa UInt64
        if !haskey(map, h)
            map[h] = Set{Hypergraph{T}}()
        end
        is_new = true
        for G2 in map[h]
            if are_isomorphic(G, G2)
                is_new = false
                break
            end
        end
        is_new && push!(map[h], G)
    end
    U = Hypergraph{T}[]
    for (_, GG) in map, G in GG
        push!(U, G)
    end
    return U
end

function _are_isomorphic(
    H1::Hypergraph{T}, H2::Hypergraph{T},
    vars1::Vector{T}, vars2::Vector{T}
) where T
    @assert length(vars1) == length(vars2)
    f = Dict{T, T}(v1 => v2 for (v1, v2) in zip(vars1, vars2))
    @assert Set{T}(f[v1] for v1 in H1.vars) == Set{T}(H2.vars)
    Set{Set{T}}(Set{T}(f[v1] for v1 in e1) for e1 in H1.edges) !=
        Set{Set{T}}(e2 for e2 in H2.edges) && return false
    return true
end

function test_isomorphism_hash1()
    H1 = Hypergraph(
        [:A, :B, :C, :D],
        [[:A, :B], [:A, :C], [:C, :D]],
    )

    H2 = Hypergraph(
        [:X, :Y, :Z, :W],
        [[:X, :Y], [:Y, :Z], [:Z, :W]],
    )

    H3 = Hypergraph(
        [:X, :Y, :Z],
        [[:X, :Y], [:Y, :Z], [:X, :Z]],
    )

    H4 = Hypergraph(
        [:X, :Y, :Z, :W],
        [[:X, :Y], [:Y, :Z], [:Z, :W], [:W, :X]],
    )

    H5 = Hypergraph(
        [:A, :B, :C, :D],
        [[:C, :B], [:D, :A], [:C, :D], [:A, :B]]
    )

    h1, _ = isomorphic_hash(H1)
    h2, _ = isomorphic_hash(H2)
    h3, _ = isomorphic_hash(H3)
    h4, _ = isomorphic_hash(H4)
    h5, _ = isomorphic_hash(H5)
    @assert h1 == h2
    @assert h1 != h3
    @assert h3 != h4
    @assert h1 != h4
    @assert h4 == h5
    @assert are_isomorphic(H1, H2)
    @assert !are_isomorphic(H1, H3)
    @assert !are_isomorphic(H1, H4)
    @assert !are_isomorphic(H3, H4)
    @assert are_isomorphic(H4, H5)
    @assert length(isomorphically_unique([H1, H2, H3, H4, H5])) == 3
end

function test_all()
    test_isomorphism_hash1()
end

end
