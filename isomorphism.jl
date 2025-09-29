module Isomorphism

using ..HypergraphWidths

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
    return Base.hash(sort(collect(values(h))))
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

    h1 = isomorphic_hash(H1)
    h2 = isomorphic_hash(H2)
    h3 = isomorphic_hash(H3)
    h4 = isomorphic_hash(H4)
    @assert h1 == h2
    @assert h1 != h3
    @assert h1 != h3
    @assert h1 != h4
end

function test_all()
    test_isomorphism_hash1()
end

end
