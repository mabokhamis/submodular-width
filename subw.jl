# import Pkg
# Pkg.add("JuMP")
# Pkg.add("Clp")
# Pkg.add("MathOptInterface")
# Pkg.add("Combinatorics")
# Pkg.add("DataStructures")
# Pkg.add("LightGraphs")

module HypergraphWidths

using JuMP, Clp, MathOptInterface
using Combinatorics
using DataStructures
using LightGraphs
using Random

"""
    Hypergraph{T}

Represents a hypergraph `H` whose vertices have type `T`. This struct has the following
fields:
    - `vars`: The set of vertices of `H`
    - `edges`: The set of hyperedges of `H`, each of which is a set of vertices
    - `tds`: The collection of tree decompositions of `H`, each of which is a collection of
    bags. Each bag in turn is a set of vertices of `H`
"""
mutable struct Hypergraph{T}
    vars::Vector{T}
    edges::Vector{Set{T}}
    tds::Vector{Vector{Set{T}}}

    # `var_index` maps a vertex `vars[i]` in `vars` to its index `i`
    var_index::Dict{T, Int}

    function Hypergraph(
        vars::Vector{T},
        edges::Vector{Vector{T}};
        tds::Vector{Vector{Set{T}}} = get_tds(edges)
    ) where T
        @assert length(unique(vars)) == length(vars) """
        Vertices of the hypergraph must be unique
        """
        var_index = Dict{T, Int}(var => i for (i, var) in enumerate(vars))
        @assert all(length(unique(edge)) == length(edge) for edge in edges) """
        Vertices of each hyperedge must be unique
        """
        edges = map(edge -> Set{T}(edge), edges)
        @assert all(reduce(union!, edges; init = Set{T}()) == Set{T}(vars)) """
        The union of all hyperedges must be equal to the set of vertices of the hypergraph
        """
        return new{T}(vars, edges, tds, var_index)
    end
end

function Base.show(io::IO, H::Hypergraph{T}) where T
    println(io, "Hypergraph with vertices: ", H.vars)
    println(io, "    and hyperedges:")
    for edge ∈ H.edges
        println(io, "        ", sort(collect(edge)))
    end
end

"""
    fractional_edge_cover(vars, edges, [verbose])

Compute the fractional edge cover number of vertices in `vars` using hyperedges in `edges`
"""
function fractional_edge_cover(
    vars::Vector{T},
    edges::Vector{Set{T}};
    verbose::Bool = false,
) where T

    # `var_edges` maps a vertex `v` in `vars` to (the indices of) hyperedges in `edges`
    # containing `v`
    var_edges = Dict{T, Vector{Int}}(v => Int[] for v in vars)
    for (i, e) in enumerate(edges), v in e
        if in(v, vars)
            push!(var_edges[v], i)
        end
    end

    # initialize a linear program
    model = Model(Clp.Optimizer)
    set_optimizer_attribute(model, "LogLevel", 0)

    n = length(vars)
    m = length(edges)

    # create a variable `λ_j` for each hyperedge `e_j` where `λ_j` represents the
    # coefficient assigned to `e_j` in a fractional edge cover of `vars`
    @variable(model, λ[1:m] >= 0.0)

    # set the objective function to be `Σ_j λ_j`
    obj = @expression(model, sum(λ[j] for j in 1:m))
    @objective(model, Min, obj)

    # for each vertex `v_i ∈ vars`, add a constraint saying that `v_i` is fractionally
    # covered by a total of at least `1.0`
    @constraint(model, con[i in 1:n], sum(λ[j] for j in var_edges[vars[i]]) >= 1.0)

    optimize!(model)

    @assert termination_status(model) == MathOptInterface.OPTIMAL

    if verbose
        println(sort(vars))
        sol = value.(λ)
        println(sol)
        println(repeat("-", 40))
    end

    obj_value = objective_value(model)
    return obj_value
end

"""
    fractional_hypertree_width(H, [verbose])

Compute the fractional hypertree width of hypergraph `H`
"""
function fractional_hypertree_width(
    H::Hypergraph{T};
    verbose::Bool = false,
) where T
    fhtw = Inf
    best_td = 0
    # for each tree decomposition `td` of `H`
    for (i, td) in enumerate(H.tds)
        # let `w` be the maximum fractional edge cover number among bags of `td`
        w = maximum(fractional_edge_cover(collect(bag), H.edges) for bag in td; init = 0.0)
        # find a `td` minimizing `w`; break ties by taking the `td` with the smallest
        # number of bags
        if w < fhtw - 1e-6 || abs(w-fhtw) <= 1e-6 && length(td) < length(H.tds[best_td])
            fhtw = w
            best_td = i
        end
    end
    if verbose
        td = H.tds[best_td]
        maximum(
            fractional_edge_cover(collect(bag), H.edges; verbose = true)
        for bag in td; init = 0.0)
    end
    return fhtw
end

"""
    zip(H, U)

Given a hypergraph `H` and a subset `U` of the vertices of `H`, encode `U` as a string of
bits. For example, `{v1, v3, v4, v8}` is encoded as the binary string `10001101`
"""
function zip(H::Hypergraph{T}, U::Set{T})::Int where T
    z = 0
    for x in U
        z |= (1 << (H.var_index[x] - 1))
    end
    return z
end

"""
    unzip(H, z)

Given a hypergraph `H` and an integer `z` representing a subset `U` of the vertices of `H`
(that was encoded using `z = zip(H, U)`), return `U`
"""
function unzip(H::Hypergraph{T}, z::Int)::Set{T} where T
    set = Set{T}()
    i = 1
    while z != 0
        if z & 1 == 1
            push!(set, H.vars[i])
        end
        i += 1
        z >>= 1
    end
    return set
end

"""
    submodular_width(H; [verbose = false])

Given a hypergraph `H` compute its submodular width. The submodular width is computed
using equation (106) in [this paper](https://arxiv.org/pdf/1612.02503v4.pdf).
"""
function submodular_width(
    H::Hypergraph{T};
    verbose::Bool = false,
) where T
    n = length(H.vars)
    N = 2 ^ n
    f(X) = sort(collect(unzip(H, X)))
    result = 0.0
    # Let `(td1, td2, ⋯, td_k)` be the (non-redundant) tree decompositions of `H`.
    # To compute the submodular width of `H`, we have to solve a linear program for each
    # combination of bags `(bag1, bag2, ⋯, bag_k)` where `bag1 ∈ td1, bag2 ∈ td2, …,`
    # `bag_k ∈ td_k` and take the maximum value across all such combinations.
    selectors = get_all_bag_selectors(H.tds)
    # selectors = Iterators.product(H.tds...,)
    println("    Final number of bag selectors: $(length(selectors))")
    counter = 0
    for β in selectors
        counter += 1
        # initialize a linear program (LP)
        model = Model(Clp.Optimizer)
        set_optimizer_attribute(model, "LogLevel", 0)

        # Let `V` be the set of vertices of `H`. For each subset `U ⊆ V`, the LP contains a
        # corresponding variable `h[U]`
        @variable(model, h[0:N-1])

        # The LP contains the constraint `h[∅] = 0`
        verbose && println("\nZero Constraint:")
        @constraint(model, h[0] == 0.0)
        verbose && println("h[$(f(0))] == 0.0")

        # For each `X ⊆ Y ⊆ V`, the LP contains a constraint `h[X] ≤ h[Y]`. These are called
        # "monotonicity constraints"
        verbose && println("\n(Elemental) Monotonicity Constraints:")
        for y = 0:n-1
            Y = N - 1
            X = Y & ~(1 << y)
            @constraint(model, h[Y] - h[X] ≥ 0.0)
            verbose && println("$(f(Y)) - $(f(X)) ≥ 0.0")
        end

        # For each `Y, Z ⊆ V` where `Y` and `Z` are not contained in one another, the LP
        # contains a constraint `h[Y] + h[Z] ≥ h[Y ∩ Z] + h[Y ∪ Z]`. These are called
        # "submodularity constraints". (Alternatively they can formulated as follows
        # using "conditional entropy" notation: `h[Y | Y ∩ Z] ≥ h[Y ∪ Z | Z]`.)

        verbose && println("\n(Elemental) Submodularity Constraints:")
        # In the submodular width case, all of the above submodularity constraints are
        # actually included in the LP. However, some of these constraints can be inferred
        # from others. Hence it suffices to include a minimal subset of the submodularity
        # constraints that is sufficient to infer all the others, which is what we do below.
        for X = 0:N-1, y = 0:n-1, z = y+1:n-1
            if (X & (1 << y) == 0) && (X & (1 << z) == 0)
                Y = X | (1 << y)
                Z = X | (1 << z)
                W = Y | (1 << z)
                @constraint(model, h[Y] + h[Z] - h[X] - h[W] ≥ 0.0)
                verbose && println("$(f(Y)) + $(f(Z)) - $(f(X)) - $(f(W)) ≥ 0.0")
            end
        end

        # For each hyperedge `e` in `H`, the LP contains a constraint `h[e] ≤ 1.0`. These
        # are called "edge-domination" constraints.
        verbose && println("\nEdge-domination Constraints:")
        for edge in H.edges
            E = zip(H, edge)
            @constraint(model, h[E] ≤ 1.0)
            verbose && println("$(f(E)) ≤ 1.0")
        end

        # The actual objective of the LP is to maximize the minimum value among
        # `h[bag1], h[bag2], …, h[bag_k]`. To that end, we add to the LP a new variable `w`
        # along with the constraints `w ≤ h[bag1], w ≤ h[bag2], …, w ≤ h[bag_k]`

        @variable(model, w >= 0.0)

        verbose && println("\nMin-target Constraints:")
        for target in β
            B = zip(H, target)
            @constraint(model, w <= h[B])
            verbose && println("w <= h[$(f(B))]")
        end

        # Finally, we set the objective of the LP to maximize `W`
        @objective(model, Max, w)

        optimize!(model)
        @assert termination_status(model) == MathOptInterface.OPTIMAL
        obj = objective_value(model)
        verbose && println("\nObjective Value: $obj")
        result = max(result, obj)
        if counter % 100 == 0
            println("        Bag selector $counter/$(length(selectors)): submodular width so far is at least $result")
        end

        verbose && println(repeat("=", 80))
    end
    return result
end

"""
    _td_from_var_order(edges, var_order)

Construct a non-redundant tree decomposition of `edges` using the variable order `var_order`
by eliminating variables in order and creating corresponding bags.
"""
function _td_from_var_order(
    edges::Set{Set{T}},
    var_order::Vector{T}
)::Set{Set{T}} where T

    edges = deepcopy(edges)
    bags = Set{Set{T}}()
    for v in var_order
        # `bag` is the union of all edges containing `v`
        bag = Set{T}()
        for e in edges
            if v in e
                union!(bag, e)
            end
        end
        isempty(bag) && continue
        # if `bag` is not empty, add it to `bags`
        push!(bags, bag)
        # `contained_edges` are all edges that are contained in `bag`
        contained_edges = Set{Set{T}}()
        for e in edges
            if issubset(e, bag)
                push!(contained_edges, e)
            end
        end
        # remove `contained_edges` from `edges`
        setdiff!(edges, contained_edges)
        # create a new edge whose variables are `bag` *minus* any variable that only
        # appears in `bag` (since these are private variables; note that `v` is one of them)
        new_edge = intersect(bag, reduce(union!, edges; init = Set{T}()))
        push!(edges, new_edge)
    end
    return bags
end

"""
    get_tds(edges)

Construct all non-redundant tree decompositions of `edges`
"""
function get_tds(edges::Vector{Vector{T}})::Vector{Vector{Set{T}}} where T
    # convert `edges` from `Vector{Vector{T}}` to `Set{Set{T}}`
    edges = Set{Set{T}}(map(edge -> Set{T}(edge), edges))
    tds = Set{Set{Set{T}}}()
    vars = collect(union(edges...,))
    # for each permutation `var_order` of the variables
    for var_order in permutations(vars)
        # construct a tree decomposition by eliminating variables using `var_order`
        td = _td_from_var_order(edges, var_order)
        # duplicate tree decompositions are automatically removed because we are storing
        # TDs in a set
        push!(tds, td)
    end
    # convert `tds` from `Set{Set{Set{T}}}` to `Vector{Vector{Set{T}}}`
    tds = map(td -> collect(td), collect(tds))
    tds = remove_subsumed_tds(tds)
    return tds
end

"""
    is_subsumed_by(td1, td2; is_td = true)

Return whether a tree decomposition `td1` is subsumed by `td2`. The optional flag `is_td`
determines whether we want to treat `td1` and `td2` as tree decompositions or as bag
selectors.
"""
function is_subsumed_by(td1::Vector{Set{T}}, td2::Vector{Set{T}}; is_td = true) where T
    return is_td ?
        all(any(issubset(bag2, bag1) for bag1 ∈ td1) for bag2 ∈ td2) :
        all(any(issubset(bag1, bag2) for bag1 ∈ td1) for bag2 ∈ td2)
end

"""
    remove_subsumed_tds(tds; is_tds = true)

Given a list of tree decompositions `tds`, remove subsumed tree decompositions and return
the resulting list. The optional flag `is_td` determines whether we want to treat `tds` as
a list of tree decompositions or as a list of bag selectors.
"""
function remove_subsumed_tds(tds::Vector{Vector{Set{T}}}; is_td = true) where T
    output_tds = Vector{Vector{Set{T}}}()
    for (i, td1) ∈ enumerate(tds)
        is_subsumed = false
        for (j, td2) ∈ enumerate(tds)
            if j != i
                if is_subsumed_by(td1, td2; is_td) && (
                        !is_subsumed_by(td2, td1; is_td) || is_subsumed_by(td2, td1; is_td) && i > j)
                    is_subsumed = true
                    break
                end
            end
        end
        if !is_subsumed
            push!(output_tds, td1)
        end
    end
    return output_tds
end

"""
    filter_selector(selector)

Given a bag selector, removed subsumed bags (i.e. that contain other bags) and return the
resulting bag selector.
"""
function filter_selector(selector::Vector{Set{T}}) where T
    new_selector = Vector{Set{T}}()
    for (i, bag1) ∈ enumerate(selector)
        is_subsumed = false
        for (j, bag2) ∈ enumerate(selector)
            if j != i
                if issubset(bag2, bag1) && (
                        bag1 != bag2 || bag1 == bag2 && i > j)
                    is_subsumed = true
                    break
                end
            end
        end
        if !is_subsumed
            push!(new_selector, bag1)
        end
    end
    return new_selector
end

"""
    get_bag_selectors(bag_selectors, td)

Given a list of `bag_selectors` and a new tree decomposition `td` that is not included in
`bag_selectors`, extend `bag_selectors` with the new `td`.
"""
function get_bag_selectors(bag_selectors::Vector{Vector{Set{T}}}, td::Vector{Set{T}}) where T
    new_selectors = Vector{Vector{Set{T}}}()
    for selector ∈ bag_selectors
        for b ∈ td
            push!(new_selectors, filter_selector([selector; b]))
        end
    end
    new_selectors = remove_subsumed_tds(new_selectors; is_td = false)
    return new_selectors
end

"""
    get_all_bag_selectors(tds)

Given a list of tree decompositions `tds`, return all possible bag selectors (not including
subsumed ones).
"""
function get_all_bag_selectors(tds::Vector{Vector{Set{T}}}) where T
    selectors = [[bag1] for bag1 ∈ first(tds)]
    println("    Number of TDs: $(length(tds))")
    println("        Creating bag selectors for TD 1/$(length(tds))")
    println("            Number of bag selectors so far: $(length(selectors))")
    for i = 2:length(tds)
        println("        Creating bag selectors for TD $i/$(length(tds))")
        selectors = get_bag_selectors(selectors, tds[i])
        println("            Number of bag selectors so far: $(length(selectors))")
    end
    return selectors
end

# println(repeat("=", 80))
# H = Hypergraph(
#     [1, 2, 3, 4],
#     [[1, 2], [2, 3], [3, 4], [4, 1]]
# )
# @show(H)
# @show(fractional_hypertree_width(H))
# @show(submodular_width(H))

# println(repeat("=", 80))
# H = Hypergraph(
#     [1, 2, 3, 4, 5],
#     [[1, 2], [2, 3], [3, 4], [4, 5], [5, 1]]
# )
# @show(H)
# @show(fractional_hypertree_width(H))
# @show(submodular_width(H))

# println(repeat("=", 80))
# H = Hypergraph(
#     [1, 2, 3, 4, 5, 6],
#     [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 1]]
# )
# @show(H)
# @show(fractional_hypertree_width(H))
# @show(submodular_width(H))

println(repeat("=", 80))
H = Hypergraph(
    ['x', 'y', 'z', 'u', 'v', 'w'],
    [
        ['x', 'w', 'z'],
        ['x', 'u', 'y'],
        ['y', 'v', 'z'],
        ['u', 'v', 'w']
    ]
)
@show(H)
@show(fractional_hypertree_width(H))
@show(submodular_width(H))

# H2 = Hypergraph(
#     [1, 2, 3, 4, 5, 6],
#     [
#         [1, 2], [1, 3], [2, 3], [4, 5], [5, 6], [4, 6],
#         [1, 5], [1, 6], [2, 4], [2, 6], [3, 4], [3, 5],
#     ]
# )
# println(fractional_hypertree_width(H2))
# println(submodular_width(H2))

# H3 = Hypergraph(
#     [1, 2, 3, 4, 5, 6, 7, 8],
#     [
#         [1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8], [8, 1],
#         [1, 5], [2, 6], [3, 7], [4, 8],
#     ],
# )

# println(submodular_width(H3))

# H4 = Hypergraph(
#     [1, 2, 3, 4, 5, 6, 7, 8],
#     [
#         [1, 3], [1, 4], [1, 5], [2, 3], [2, 4], [2, 5], [3, 6], [3, 7], [4, 6], [4, 8],
#         [5, 7], [5, 8]
#     ],
# )
# println(fractional_hypertree_width(H4))
# println(submodular_width(H4))

# H5 = Hypergraph(
#     [1, 2, 3, 4, 5, 6, 7, 8 ,9, 10, 11, 12],
#     [
#         [1, 2], [2, 4], [1, 3], [3, 4],
#         [4, 5], [5, 7], [4, 6], [6, 7],
#         [7, 9], [9, 10], [7, 8], [8, 10],
#         [10, 11], [11, 1], [10, 12], [12, 1]
#     ]
# )

# println(fractional_hypertree_width(H5))
# println(submodular_width(H5))

# α = 8/13

# Ha = Hypergraph(
#     [1, 2, 3, 4],
#     [[1, 2], [2, 3], [3, 4], [4, 1]];
#     weights = [1+α, 1+α, 1+α, 1+α]
# )

# println(submodular_width(Ha))

# Hb = Hypergraph(
#     [1, 2, 3, 4, 5, 6, 7, 8],
#     [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 7], [7, 8], [8, 1]];
#     weights = [2-α, 2-α, 2-α, 2-α, 2-α, 2-α, 2-α, 2-α]
# )
# println(submodular_width(Hb))

# Hc = Hypergraph(
#     [1, 2, 3, 4, 5],
#     [[1, 2], [2, 3], [3, 4], [4, 5], [5, 1]];
#     weights = [1+α, 1+α, 1+α, 2-α, 2-α]
# )

# println(submodular_width(Hc))

# Hd = Hypergraph(
#     [1, 2, 3, 4, 5, 6],
#     [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 1]];
#     weights = [1+α, 1+α, 2-α, 2-α, 2-α, 2-α]
# )

# println(submodular_width(Hd))

# tds = [
#     Set{Int64}[Set([5, 6, 7, 2]), Set([5, 7, 8, 3, 1]), Set([4, 8, 3, 1]), Set([5, 7, 2, 3, 1])],
#     Set{Int64}[Set([6, 2, 3, 1]), Set([5, 4, 7, 3, 1]), Set([5, 4, 7, 8]), Set([5, 6, 7, 3, 1])],
#     # Set{Int64}[Set([4, 8, 3, 1]), Set([5, 6, 2, 8, 3]), Set([5, 2, 8, 3, 1]), Set([6, 7, 8, 3])],
#     Set{Int64}[Set([6, 8, 3, 1]), Set([6, 2, 3, 1]), Set([4, 8, 3, 1]), Set([5, 6, 8, 1]), Set([6, 7, 8, 3])],
#     Set{Int64}[Set([4, 7, 2, 3]), Set([5, 6, 7, 2]), Set([5, 4, 7, 2]), Set([5, 4, 2, 1]), Set([5, 4, 7, 8])],
#     # Set{Int64}[Set([4, 6, 2, 8, 3]), Set([5, 4, 2, 1]), Set([5, 4, 6, 2, 8]), Set([6, 7, 8, 3])],
# ]

# pretty_tds = [map(bag -> collect(bag), td) for td in tds]
# for td in pretty_tds
#     for bag in td
#         sort!(bag)
#     end
#     sort!(td)
#     println(td)
# end
# println()

# HC = Hypergraph(
#     [1, 2, 3, 4, 5, 6, 7, 8],
#     [
#         [1, 2], [2, 3], [3, 4], [4, 1],
#         [5, 6], [6, 7], [7, 8], [8, 5],
#         [1, 5], [2, 6], [3, 7], [4, 8]
#     ];
#     tds = tds
# )

# # println(fractional_hypertree_width(HC))
# println(submodular_width(HC))

function get_multivariate_extension(H::Hypergraph{T}, extras::Vector{T}) where T
    @assert isempty(extras ∩ H.vars) && length(extras) == length(H.edges)
    m = length(H.edges)
    E = Hypergraph{T}[]
    for p in permutations(1:m)
        # @show(p)
        edges = deepcopy(H.edges)
        for (i, k) in enumerate(p)
            union!(edges[k], extras[2:min(i, m-1)])
        end
        edges = map(edge -> collect(edge), edges)
        push!(E, Hypergraph(H.vars ∪ extras[2:end-1], edges))
    end
    return E
end

function fractional_hypertree_width(E::Vector{Hypergraph{T}}) where T
    m = 0.0
    witness_reported = false
    for (i, H) in enumerate(E)
        w = fractional_hypertree_width(H)
        if abs(w - 2.0) < 1e-6 && !witness_reported
            @warn "$H"
            witness_reported = true
        end
        m = max(m, w)
        # println("$i: $m")
    end
    return m
end

# ----------------------------------------------------------------------------

# H = Hypergraph(
#     [1, 2, 3, 4],
#     [[1, 2, 3], [2, 3, 4], [3, 4, 1], [4, 1, 2]],
# )

# E = get_multivariate_extension(H, [-1, -2, -3, -4])
# println(fractional_hypertree_width(H))
# println(fractional_hypertree_width(E))


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



# ======================================================================================== #
# Ahmet's queries:
# ================

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

end
