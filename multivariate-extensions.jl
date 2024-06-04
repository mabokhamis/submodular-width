"""
    MultivariateExtensions

This module contains tools for computing the fractional hypertree width of a _multivariate
# extension_ of a query. Multivariate extensions are defined in [this
# paper](https://arxiv.org/abs/2312.09331)
"""
module MultivariateExtensions

using Combinatorics

using ..HypergraphWidths

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
        push!(E, Hypergraph(H.vars ∪ Z[2:end-1], edges))
    end
    return E
end

function HypergraphWidths.fractional_hypertree_width(E::Vector{Hypergraph{T}}) where T
    m = 0.0
    witness_reported = false
    println("    Computing fractional hypertree Width of multivariate extensions:")
    for (i, H) in enumerate(E)
        w = fractional_hypertree_width(H)
        if abs(w - 2.0) < 1e-6 && !witness_reported
            @warn "$H"
            witness_reported = true
        end
        m = max(m, w)
        println("        Extension $i/$(length(E)): FHTW so far: $m")
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


# # 3-path
# H = Hypergraph(
#     [:A, :B, :C, :D],
#     [[:A], [:A, :B], [:B, :C], [:C, :D], [:D]],
# )

# E = get_multivariate_extension(H, [:Z1, :Z2, :Z3, :Z4, :Z5])
# @warn "$(length(E))"
# println("RESULT:", fractional_hypertree_width(E))

# 2-path
# H = Hypergraph(
#     [:A, :B, :C],
#     [[:A], [:A, :B], [:B, :C], [:C]],
# )

# E = get_multivariate_extension(H, [:Z1, :Z2, :Z3, :Z4])
# @warn "$(length(E))"
# println("RESULT:", fractional_hypertree_width(E))

# # 1-path
# H = Hypergraph(
#     [:A, :B],
#     [[:A], [:A, :B], [:B]],
# )

# E = get_multivariate_extension(H, [:Z1, :Z2, :Z3])
# @warn "$(length(E))"
# println("RESULT:", fractional_hypertree_width(E))



# ------------------------------------------------------------------------------------
# Ahmet's queries:
# ================

# # triangle + one atom
# H = Hypergraph(
#     [:A, :B, :C],
#     [[:A, :B], [:B, :C], [:C, :A], [:C]],
# )

# E = get_multivariate_extension(H, [:Z1, :Z2, :Z3, :Z4])
# @warn "$(length(E))"
# println("RESULT:", fractional_hypertree_width(E))
# # Result: 5/3

# # ----------------------------------------------------------------------------

# # Q_1
# H = Hypergraph(
#     [:A, :B, :C, :D, :E, :F],
#     [[:A, :B], [:B, :C], [:C, :D], [:B, :E], [:C, :F]],
# )

# E = get_multivariate_extension(H, [:Z1, :Z2, :Z3, :Z4, :Z5])
# @warn "$(length(E))"
# println("RESULT:", fractional_hypertree_width(E))
# # Result: 1.5

# # ----------------------------------------------------------------------------

# # Q_1'
# H = Hypergraph(
#     [:A, :B, :C, :D, :E, :E2, :F],
#     [[:A, :B], [:B, :C], [:C, :D], [:B, :E], [:B, :E2], [:C, :F]],
# )

# E = get_multivariate_extension(H, [:Z1, :Z2, :Z3, :Z4, :Z5, :Z6])
# @warn "$(length(E))"
# println("RESULT:", fractional_hypertree_width(E))
# # Result: 1.5

# # ----------------------------------------------------------------------------

# # Q_2
# H = Hypergraph(
#     [:A, :B, :C],
#     [[:A, :B, :C], [:A, :B], [:B, :C]],
# )

# E = get_multivariate_extension(H, [:Z1, :Z2, :Z3])
# @warn "$(length(E))"
# println("RESULT:", fractional_hypertree_width(E))
# # Result: 1.5

# # ----------------------------------------------------------------------------

# # Q_2'
# H = Hypergraph(
#     [:A, :B, :C, :D],
#     [[:A, :B, :C, :D], [:A, :B], [:B, :C], [:C, :D]],
# )

# E = get_multivariate_extension(H, [:Z1, :Z2, :Z3, :Z4])
# @warn "$(length(E))"
# println("RESULT:", fractional_hypertree_width(E))
# # Result: 1.5

# # ----------------------------------------------------------------------------

# # Q_2''
# H = Hypergraph(
#     [:A, :B, :C, :D],
#     [[:A, :B, :C, :D], [:A, :B], [:B, :C], [:B, :D]],
# )

# E = get_multivariate_extension(H, [:Z1, :Z2, :Z3, :Z4])
# @warn "$(length(E))"
# println("RESULT:", fractional_hypertree_width(E))
# # Result: 5/3

# # ----------------------------------------------------------------------------

# # Q_3
# H = Hypergraph(
#     [:A, :B, :C],
#     [[:A, :B, :C], [:A], [:B], [:C]],
# )

# E = get_multivariate_extension(H, [:Z1, :Z2, :Z3, :Z4])
# @warn "$(length(E))"
# println("RESULT:", fractional_hypertree_width(E))
# # Result: 5/3

# # ----------------------------------------------------------------------------

# # Q_3'
# H = Hypergraph(
#     [:A, :B, :C, :D],
#     [[:A, :B, :C, :D], [:A], [:B], [:C], [:D]],
# )

# E = get_multivariate_extension(H, [:Z1, :Z2, :Z3, :Z4, :Z5])
# @warn "$(length(E))"
# println("RESULT:", fractional_hypertree_width(E))
# # Result: 1:75

#=========================================================================================#

end
