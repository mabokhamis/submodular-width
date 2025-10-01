"""
    MultivariateExtensions

This module contains tools for computing the fractional hypertree width of a _multivariate
# extension_ of a query. Multivariate extensions are defined in [this
# paper](https://arxiv.org/abs/2312.09331)
"""
module MultivariateExtensions

using Combinatorics

using ..HypergraphWidths
using ..Isomorphism

"""
    get_multivariate_extension(H, Z)

Compute the multivariate extensions of a given hypergraph `H`. `Z` is a vector of new
variables to be used for the multivariate extension. See [this #
    paper](https://arxiv.org/abs/2312.09331) for more details.
"""
function get_multivariate_extension(H::Hypergraph{T}, Z::Vector{T}) where T
    @assert isempty(Z ∩ H.vars) && length(Z) == length(H.edges) """
    The extension variables `Z` must be disjoint from the variables of the hypergraph `H`.
    The number of extension variables must be equal to the number of edges in the hypergraph.
    """
    m = length(H.edges)
    E = Hypergraph{T}[]
    println("    Computing multivariate extensions:")
    for p in permutations(1:m)
        println("        $p")
        edges = deepcopy(H.edges)
        for (i, k) in enumerate(p)
            union!(edges[k], Z[2:min(i, m-1)])
        end
        edges = map(edge -> collect(edge), edges)
        G = Hypergraph(H.vars ∪ Z[2:end-1], edges)
        push!(E, simplify_hypergraph(G))
    end
    return isomorphically_unique(E)
end

function HypergraphWidths.fractional_hypertree_width(E::Vector{Hypergraph{T}}) where T
    m = 0.0
    println("    Computing fractional hypertree Width of multivariate extensions:")
    for (i, H) in enumerate(E)
        w = fractional_hypertree_width(H)
        if w < m - 1e-6
            @warn "New Hypergraph found:\n$H\nwith width $w < $m"
        end
        m = max(m, w)
        println("        Extension $i/$(length(E)): FHTW so far: $m")
    end
    return m
end

function HypergraphWidths.submodular_width(E::Vector{Hypergraph{T}}) where T
    m = 0.0
    println("    Computing submodular width of multivariate extensions:")
    for (i, H) in enumerate(E)
        w = submodular_width(H)
        if w < m - 1e-6
            @warn "New Hypergraph found:\n$H\nwith width $w < $m"
        end
        m = max(m, w)
        println("        Extension $i/$(length(E)): SUBW so far: $m")
    end
    return m
end

#==========================================================================================#
# Testcases:
# ----------

function test_triangle()
    println(repeat("=", 80))
    H = Hypergraph(
        [:A, :B, :C],
        [[:A, :B], [:B, :C], [:C, :A]],
    )
    println(H)
    E = get_multivariate_extension(H, [:Z1, :Z2, :Z3])
    @assert fractional_hypertree_width(H) ≈ 1.5
    @assert fractional_hypertree_width(E) ≈ 1.5
end

function test_LW4()
    println(repeat("=", 80))
    H = Hypergraph(
        [1, 2, 3, 4],
        [[1, 2, 3], [2, 3, 4], [3, 4, 1], [4, 1, 2]],
    )
    println(H)
    E = get_multivariate_extension(H, [-1, -2, -3, -4])
    @assert fractional_hypertree_width(H) ≈ 4/3
    @assert fractional_hypertree_width(E) ≈ 1.5
end

function test_2path()
    println(repeat("=", 80))
    H = Hypergraph(
        [:A, :B, :C],
        [[:A, :B], [:B, :C]],
    )
    println(H)
    E = get_multivariate_extension(H, [:Z1, :Z2])
    @assert fractional_hypertree_width(H) ≈ 1
    @assert fractional_hypertree_width(E) ≈ 1
end

function test_bowtie()
    println(repeat("=", 80))
    H = Hypergraph(
        [:A, :B],
        [[:A], [:A, :B], [:B]],
    )
    println(H)
    E = get_multivariate_extension(H, [:Z1, :Z2, :Z3])
    @assert fractional_hypertree_width(H) ≈ 1
    @assert fractional_hypertree_width(E) ≈ 1.5
end

function test_2path_with_endpoints()
    println(repeat("=", 80))
    H = Hypergraph(
        [:A, :B, :C],
        [[:A], [:A, :B], [:B, :C], [:C]],
    )
    println(H)
    E = get_multivariate_extension(H, [:Z1, :Z2, :Z3, :Z4])
    @assert fractional_hypertree_width(H) ≈ 1
    @assert fractional_hypertree_width(E) ≈ 2
end

function test_all()
    test_triangle()
    test_LW4()
    test_2path()
    test_bowtie()
    test_2path_with_endpoints()
end

H = Hypergraph(
    [:A, :B, :C, :D, :E, :F],
    [[:A, :B], [:B, :C], [:C, :D], [:D, :E], [:E, :F]],
)

query_name = "5path"

io = open("$(query_name).jl", "w")
println(io, "using ..HypergraphWidths")
println(io, "io = open(\"$(query_name).log\", \"a\")")
println(io, "println(io, repeat(\"#\", 80))")

E = get_multivariate_extension(H, [:Z1, :Z2, :Z3, :Z4, :Z5])
for (i, H) in enumerate(E)
    println(io, repeat("#", 80))
    println(io, "# Extension $i/$(length(E)):")
    println(io, "H = $(H)")
    println(io, "print(io, \"$i/$(length(E)): \t\")\nprintln(io, submodular_width(H))")
    println(io, "flush(io)")
end
println(io, "close(io)")
close(io)

end
