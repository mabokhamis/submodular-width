"""
    HypergraphWidths

This module contains tools for computing the fractional hypertree width and submodular width
of a query. In addition, it supports the FD-aware and degree-aware variants of those width
measures.
"""
module HypergraphWidths

using JuMP
using Clp
using MathOptInterface
using Combinatorics
using DataStructures

export Hypergraph, FD, fractional_edge_cover, fractional_hypertree_width, submodular_width,
    get_tds, get_trivial_tds

"""
    Hypergraph{T}

Represents a hypergraph `H` whose vertices have type `T`. This struct has the following
fields:
    - `vars`: The set of vertices of `H`
    - `edges`: The set of hyperedges of `H`, each of which is a set of vertices
    - `weights`: The weights of the hyperedges of `H`. By default, all weights are `1.0`.
      When computing bounds with proper degree constraints (like the polymatroid bound or
      the degree-aware submodular width), the weight of a hyperedge represents the log of
      the size of the corresponding relation. (The weight can be ∞ if the relation is
      infinite, e.g. `x + y = z`)
    - `tds`: The collection of tree decompositions of `H`, each of which is a collection of
      bags. Each bag in turn is a set of vertices of `H`
"""
mutable struct Hypergraph{T}
    vars::Vector{T}
    edges::Vector{Set{T}}
    weights::Vector{Float64}
    tds::Vector{Vector{Set{T}}}

    # `_var_index` maps a vertex `vars[i]` in `vars` to its index `i`
    _var_index::Dict{T, Int}

    # `_var_edges` maps a vertex `v` in `vars` to (the indices of) hyperedges in `edges`
    # containing `v`
    _var_edges::Dict{T, Set{Int}}

    """
        Hypergraph(vars, edges; weights = ones(length(edges)), tds = get_tds(edges))

    Construct a hypergraph `H` with vertices `vars` and hyperedges `edges`. Optional
    `weights` and `tds` can be provided. By default, the weights of the hyperedges are all
    `1.0`, and the vector oftree decompositions `tds` is empty.
    """
    function Hypergraph(
        vars::Vector{T},
        edges::Union{Vector{Vector{T}},Vector{Set{T}}};
        weights::Vector{Float64} = ones(length(edges)),
        tds::Vector{Vector{Set{T}}} = Vector{Vector{Set{T}}}()
        # Alternatively, use `tds = get_trivial_tds(edges)` to skip the above expensive
        # computation of TDs, if TDs are not needed (e.g. for the polymatroid bound)
    ) where T
        @assert length(unique(vars)) == length(vars) """
        Vertices of the hypergraph must be unique
        """
        @assert all(length(unique(edge)) == length(edge) for edge in edges) """
        Vertices of each hyperedge must be unique
        """
        edges = map(edge -> Set{T}(edge), edges)
        @assert all(reduce(union!, edges; init = Set{T}()) == Set{T}(vars)) """
        The union of all hyperedges must be equal to the set of vertices of the hypergraph
        """
        @assert length(weights) == length(edges) """
        The number of weights must be equal to the number of hyperedges
        """
        @assert all(w ≥ 0.0 for w in weights) """
        Weights must be non-negative
        """

        _var_index = Dict{T, Int}(var => i for (i, var) in enumerate(vars))

        _var_edges = Dict{T, Set{Int}}(v => Set{Int}() for v in vars)
        for (i, e) in enumerate(edges), v in e
            push!(_var_edges[v], i)
        end
        return new{T}(vars, edges, weights, tds, _var_index, _var_edges)
    end
end

# Copy a hypergraph
function Base.copy(H::Hypergraph{T}) where T
    if length(H.tds) == 0
        tds = Vector{Vector{Set{T}}}()
    else
        tds = [[copy(bag) for bag in td] for td in H.tds]
    end
    return Hypergraph(
        copy(H.vars),
        [copy(E) for E in H.edges];
        weights = copy(H.weights),
        tds
    )
end

function Base.show(io::IO, H::Hypergraph{T}) where T
    println(io, "Hypergraph with vertices: ", H.vars)
    println(io, "    and hyperedges:")
    for (i, edge) ∈ enumerate(H.edges)
        println(io, "        ", sort(collect(edge)), " with weight $(H.weights[i])")
    end
end

# if the TDs are missing in the hypergraph, initize them to contain all possible TDs
function initialize_tds_if_missing(H)
    if length(H.tds) == 0
        H.tds = get_tds(H.edges)
    end
end

"""
    FD{T}

A functional dependency (FD) `X → Y` where `X` and `Y` are sets of vertices of the query's
hypergraph
"""
struct FD{T}
    X::Set{T}
    Y::Set{T}

    function FD(X::Union{Vector{T},Set{T}}, Y::Union{Vector{T},Set{T}}) where T
        @assert length(unique(X)) == length(X) """
        In an FD `$X → $Y`, the variables in `$X` must be unique
        """
        @assert length(unique(Y)) == length(Y) """
        In an FD `$X → $Y`, the variables in `$Y` must be unique
        """
        @assert isdisjoint(X, Y) """
        In an FD `$X → $Y`, the sets `$X and `$Y` must be disjoint
        """
        return new{T}(Set{T}(X), Set{T}(Y) ∪ Set{T}(X))
    end
end

function Base.show(io::IO, fd::FD{T}) where T
    X = sort(collect(fd.X))
    Y = sort(collect(setdiff(fd.Y, fd.X)))
    print(io, "$X → $Y")
end

function Base.show(io::IO, fds::Vector{FD{T}}) where T
    println(io, "Functional Dependencies:")
    for fd ∈ fds
        println(io, "    $fd")
    end
end

"""
    DC{T}

A Degree Constraint DC is represented by a triple `(X, Y, n)` where `X` and `Y` are sets of
variables are n is a non-negative number. It means:
```
log max deg(Y|X) ≤ n
```
NOTE that `n` is on log-scale: it is an upper bound on the **log** of the maximum
degree of `Y` given `X`
"""
struct DC{T}
    X::Set{T}
    Y::Set{T}
    n::Float64

    function DC(X::Union{Vector{T},Set{T}}, Y::Union{Vector{T},Set{T}}, n::Number) where T
        @assert length(unique(X)) == length(X) """
        In a DC `($X, $Y, $n)`, the variables in `$X` must be unique
        """
        @assert length(unique(Y)) == length(Y) """
        In a DC `($X, $Y, $n)`, the variables in `$Y` must be unique
        """
        @assert isdisjoint(X, Y) """
        In a DC `($X, $Y, $n)`, the sets `$X and `$Y` must be disjoint
        """
        @assert n ≥ 0.0 """
        In a DC `($X, $Y, $n)`, the number `$n` must be non-negative
        """
        return new{T}(Set{T}(X), Set{T}(Y) ∪ Set{T}(X), Float64(n))
    end

    # An FD is a special case of a DC where n = 0
    function DC(fd::FD{T}) where T
        return new{T}(fd.X, fd.Y, 0.0)
    end
end

function Base.show(io::IO, dc::DC{T}) where T
    X = sort(collect(dc.X))
    Y = sort(collect(setdiff(dc.Y, dc.X)))
    print(io, "log max deg($Y | $X) ≤ $(dc.n)")
end

function Base.show(io::IO, dcs::Vector{DC{T}}) where T
    println(io, "Degree Constraints:")
    for dc ∈ dcs
        println(io, "    $dc")
    end
end

"""
    fractional_edge_cover(H, target_vars = H.vars; [verbose])

Compute the fractional edge cover number of a given set of vertices `target_vars` in a
target hypergraph `H`
"""
function fractional_edge_cover(
    H::Hypergraph{T},
    target_vars::Vector{T} = H.vars;
    verbose::Bool = false,
) where T
    @assert target_vars ⊆ H.vars """
    `fractional_edge_cover(H, target_vars)` expects `target_vars` to be a subset of the
    vertices of the given hypergraph `H`
    """

    # initialize a linear program
    model = Model(Clp.Optimizer)
    set_optimizer_attribute(model, "LogLevel", 0)

    n = length(target_vars) # number of constraints
    m = length(H.edges)     # number of variables

    # create a variable `λ_j` for each hyperedge `e_j` where `λ_j` represents the
    # coefficient assigned to `e_j` in a fractional edge cover of `vars`
    @variable(model, λ[1:m] >= 0.0)

    # set the objective function to be `Σ_j weight_j * λ_j`
    obj = @expression(model, sum(H.weights[j] * λ[j] for j in 1:m))
    @objective(model, Min, obj)

    # for each vertex `v_i ∈ vars`, add a constraint saying that `v_i` is fractionally
    # covered by a total of at least `1.0`
    @constraint(model, con[i in 1:n], sum(λ[j] for j in H._var_edges[target_vars[i]]) >= 1.0)

    optimize!(model)

    @assert termination_status(model) == MathOptInterface.OPTIMAL

    if verbose
        sol = value.(λ)
        println(sol)
        println(repeat("-", 40))
    end

    obj_value = objective_value(model)
    return obj_value
end

"""
    fractional_hypertree_width(H, [verbose])

Compute the fractional hypertree width of hypergraph `H`. If the vector of
tree decompositions of `H` is empty, it is updated to contain all tree decompositions
of `H`.
"""
function fractional_hypertree_width(
    H::Hypergraph{T};
    verbose::Bool = false,
) where T
    initialize_tds_if_missing(H)
    fhtw = Inf
    best_td = 0
    # for each tree decomposition `td` of `H`
    for (i, td) in enumerate(H.tds)
        # let `w` be the maximum fractional edge cover number among bags of `td`
        w = maximum(fractional_edge_cover(H, collect(bag)) for bag in td; init = 0.0)
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
            fractional_edge_cover(H, collect(bag); verbose = true)
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
        z |= (1 << (H._var_index[x] - 1))
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
    submodular_width(H; fds = FD[], dcs = DC[], verbose = false)

Given a hypergraph `H` compute its submodular width. The submodular width is computed
using equation (106) in [this paper](https://arxiv.org/pdf/1612.02503v4.pdf).

 - `fds` is an optional list of FDs
 - `dcs` is an optional list of DCs

 If the vector of tree decompositions is empty, it is updated to contain all possible
 tree decompositions of the hypergraph.
"""
function submodular_width(
    H::Hypergraph{T};
    fds::Vector{FD{T}} = FD{T}[],
    dcs::Vector{DC{T}} = DC{T}[],
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
    initialize_tds_if_missing(H)
    selectors = _get_bag_selectors(H.tds)
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
        for (i, edge) in enumerate(H.edges)
            isinf(H.weights[i]) && continue
            E = zip(H, edge)
            @constraint(model, h[E] ≤ H.weights[i])
            verbose && println("$(f(E)) ≤ $(H.weights[i])")
        end

        # Convert FDs into DCs
        dcs = copy(dcs)
        append!(dcs, DC(fd) for fd in fds)

        # For each degree constraint `log max deg(Y|X) ≤ n` in `dcs`, the LP contains a
        # constraint:
        # h[Y] - h[X] ≤ n
        # As a special case, for each FD `X → Y`, the LP contains a constraint:
        # h[Y] - h[X] ≤ 0
        verbose && println("\nDC Constraints:")
        for dc ∈ dcs
            @assert any(dc.Y ⊆ E for E in H.edges) """
            DC variables must be a contained in a hyperedge of the hypergraph. The following
            DC does not satisfy this condition: $dc
            """
            isinf(dc.n) && continue
            X = zip(H, dc.X)
            Y = zip(H, dc.Y)
            @constraint(model, h[Y] - h[X] ≤ dc.n)
            verbose && println("$(f(Y)) - $(f(X)) ≤ $(dc.n)")
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
        verbose && println("\nObjective: Maximize w")

        optimize!(model)
        @assert termination_status(model) == MathOptInterface.OPTIMAL
        obj = objective_value(model)
        verbose && println("\nOptimal Objective Value: $obj")
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

Construct all non-redundant tree decompositions of `edges`. Remove tree decompositions that
are "subsumed" by others. If a tree decomposition `td1` is subsumed by `td2`, then including
`td1` in the computation of fractional hypertree width or submodular width is redundant.
"""
function get_tds(edges::Union{Vector{Vector{T}},Vector{Set{T}}})::Vector{Vector{Set{T}}} where T
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
    tds = _remove_subsumed_tds(tds)
    return tds
end

"""
    get_trivial_tds(edges)
    get_trivial_tds(vars)

Given a list of edges (or vertices) of a hypergraph, return a list of tree decompositions
that contain only the trivial tree decomposition consisting of a single bag containing all
vertices of the hypergraph
"""
function get_trivial_tds(
    edges::Union{Vector{Vector{T}},Vector{Set{T}}}
)::Vector{Vector{Set{T}}} where T
    vars = reduce(union!, edges; init = Set{T}())
    return get_trivial_tds(vars)
end
function get_trivial_tds(
    vars::Union{Vector{T},Set{T}}
)::Vector{Vector{Set{T}}} where T
    return [[Set{T}(vars)]]
end

"""
    _is_subsumed_by(td1, td2; is_td = true)

Return whether a tree decomposition `td1` is subsumed by `td2`. The optional flag `is_td`
determines whether we want to treat `td1` and `td2` as tree decompositions or as bag
selectors.
"""
function _is_subsumed_by(td1::Vector{Set{T}}, td2::Vector{Set{T}}; is_td = true) where T
    return is_td ?
        all(any(issubset(bag2, bag1) for bag1 ∈ td1) for bag2 ∈ td2) :
        all(any(issubset(bag1, bag2) for bag1 ∈ td1) for bag2 ∈ td2)
end

"""
    _remove_subsumed_tds(tds; is_tds = true)

Given a list of tree decompositions `tds`, remove subsumed tree decompositions and return
the resulting list. The optional flag `is_td` determines whether we want to treat `tds` as
a list of tree decompositions or as a list of bag selectors.
"""
function _remove_subsumed_tds(tds::Vector{Vector{Set{T}}}; is_td = true) where T
    output_tds = Vector{Vector{Set{T}}}()
    for (i, td1) ∈ enumerate(tds)
        is_subsumed = false
        for (j, td2) ∈ enumerate(tds)
            if j != i
                if _is_subsumed_by(td1, td2; is_td) && (
                        !_is_subsumed_by(td2, td1; is_td) || _is_subsumed_by(td2, td1; is_td) && i > j)
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
    _filter_selector(selector)

Given a bag selector, removed subsumed bags (i.e. that contain other bags) and return the
resulting bag selector.
"""
function _filter_selector(selector::Vector{Set{T}}) where T
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
    _extend_bag_selectors(bag_selectors, td)

Given a list of `bag_selectors` and a new tree decomposition `td` that is not included in
`bag_selectors`, extend `bag_selectors` with the new `td`.
"""
function _extend_bag_selectors(bag_selectors::Vector{Vector{Set{T}}}, td::Vector{Set{T}}) where T
    new_selectors = Vector{Vector{Set{T}}}()
    for selector ∈ bag_selectors
        for b ∈ td
            push!(new_selectors, _filter_selector([selector; b]))
        end
    end
    new_selectors = _remove_subsumed_tds(new_selectors; is_td = false)
    return new_selectors
end

"""
    _get_bag_selectors(tds)

Given a list of tree decompositions `tds`, return all possible bag selectors (not including
subsumed ones).
"""
function _get_bag_selectors(tds::Vector{Vector{Set{T}}}) where T
    selectors = [[bag1] for bag1 ∈ first(tds)]
    println("    Number of TDs: $(length(tds))")
    println("        Creating bag selectors for TD 1/$(length(tds))")
    println("            Number of bag selectors so far: $(length(selectors))")
    for i = 2:length(tds)
        println("        Creating bag selectors for TD $i/$(length(tds))")
        selectors = _extend_bag_selectors(selectors, tds[i])
        println("            Number of bag selectors so far: $(length(selectors))")
    end
    return selectors
end

"""
    polymatroid_bound(H; fds = FD[], dcs = DC[], verbose = false)

Compute the polymatroid bound of a hypergraph `H`. The polymatroid bound is a special case
of the submodular width where we only use a single tree decomposition with a single bag
containing all variables of the hypergraph.

- `fds` is an optional list of FDs
- `dcs` is an optional list of DCs

*NOTE* The default weights of the hyperedges are all `1.0`. However, when computing the
polymatroid bound with proper DC constraints, those weights should be set to the log of the
sizes of the corresponding relations. (The weight can be ∞ if the relation is infinite, e.g.
`x + y = z`)
"""
function polymatroid_bound(
    H::Hypergraph{T};
    fds::Vector{FD{T}} = FD{T}[],
    dcs::Vector{DC{T}} = DC{T}[],
    verbose::Bool = false,
) where T
    # The polymatroid bound is a special case of the submodular width where we only use a
    # single tree decomposition with a single bag containing all variables of the hypergraph
    H = copy(H)
    H.tds = get_trivial_tds(H.vars)
    return submodular_width(H; fds, dcs, verbose)
end

#==========================================================================================#
# Testcases:
# ----------

# 4-cycle query
function test_4cycle()
    println(repeat("=", 80))
    H = Hypergraph(
        [1, 2, 3, 4],
        [[1, 2], [2, 3], [3, 4], [4, 1]]
    )
    @show(H)
    fhtw = fractional_hypertree_width(H)
    @show(fhtw)
    @assert fhtw ≈ 2.0
    subw = submodular_width(H)
    @show(subw)
    @assert subw ≈ 1.5
end

# 4-cycle with FDs
function test_4cycle_with_fds()
    println(repeat("=", 80))
    H = Hypergraph(
        [1, 2, 3, 4],
        [[1, 2], [2, 3], [3, 4], [4, 1]]
    )
    @show(H)
    fhtw = fractional_hypertree_width(H)
    @show(fhtw)
    @assert fhtw ≈ 2.0
    fds = FD{Int}[
        FD([1], [2]),
        FD([3], [2]),
    ]
    println(fds)
    subw = submodular_width(H; fds)
    @show(subw)
    @assert subw ≈ 1.0
end

# 5-cycle:
function test_5cycle()
    println(repeat("=", 80))
    H = Hypergraph(
        [1, 2, 3, 4, 5],
        [[1, 2], [2, 3], [3, 4], [4, 5], [5, 1]]
    )
    @show(H)
    fhtw = fractional_hypertree_width(H)
    @show(fhtw)
    @assert fhtw ≈ 2.0
    subw = submodular_width(H)
    @show(subw)
    @assert subw ≈ 5/3
end

# 5-cycle with FDs
function test_5cycle_with_fds()
    println(repeat("=", 80))
    H = Hypergraph(
        [1, 2, 3, 4, 5],
        [[1, 2], [2, 3], [3, 4], [4, 5], [5, 1]]
    )
    @show(H)
    fhtw = fractional_hypertree_width(H)
    @show(fhtw)
    @assert fhtw ≈ 2.0
    fds = FD{Int}[
        FD([1], [5]),
        FD([5], [1]),
    ]
    println(fds)
    subw = submodular_width(H; fds)
    @show(subw)
    @assert subw ≈ 1.5
end

# 6-cycle
function test_6cycle()
    println(repeat("=", 80))
    H = Hypergraph(
        [1, 2, 3, 4, 5, 6],
        [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 1]]
    )
    @show(H)
    fhtw = fractional_hypertree_width(H)
    @show(fhtw)
    @assert fhtw ≈ 2.0
    subw = submodular_width(H)
    @show(subw)
    @assert subw ≈ 5/3
end

# 6-cycle with FDs:
function test_6cycle_with_fds()
    println(repeat("=", 80))
    H = Hypergraph(
        [1, 2, 3, 4, 5, 6],
        [[1, 2], [2, 3], [3, 4], [4, 5], [5, 6], [6, 1]]
    )
    @show(H)
    fhtw = fractional_hypertree_width(H)
    @show(fhtw)
    @assert fhtw ≈ 2.0
    fds = FD{Int}[
        FD([2], [3]),
        FD([4], [5]),
        FD([6], [1]),
    ]
    println(fds)
    subw = submodular_width(H; fds)
    @show(subw)
    @assert subw ≈ 1.5
    end

# Example 6 on page 28 here: https://arxiv.org/pdf/1712.07880
function test_example_6()
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
    fhtw = fractional_hypertree_width(H)
    @show(fhtw)
    @assert fhtw ≈ 2.0
    fds = FD{Char}[
        FD(['x', 'y'], ['u']),
        FD(['y', 'u'], ['x']),
        FD(['u', 'x'], ['y']),

        FD(['z', 'y'], ['v']),
        FD(['y', 'v'], ['z']),
        FD(['v', 'z'], ['y']),

        FD(['x', 'z'], ['w']),
        FD(['z', 'w'], ['x']),
        FD(['w', 'x'], ['z']),
    ]
    println(fds)
    subw_no_fds = submodular_width(H)
    println("Submodular width *WITHOUT* FDs: $subw_no_fds\n")      # 1.75
    @assert subw_no_fds ≈ 1.75
    subw_fds = submodular_width(H; fds)
    println("Submodular width *WITH*    FDs: $subw_fds\n")         # 1.5
    @assert subw_fds ≈ 1.5
end

# Example from Figure 1 on page 8 here: https://arxiv.org/pdf/1604.00111
# Q :- R(x, y), S(y, z), T(z, u), xz --> u, yu --> x
function test_polymatroid_bound1()
    println(repeat("=", 80))
    H = Hypergraph(
        ['x', 'y', 'z', 'u'],
        [
            ['x', 'y'],      # R(x, y)
            ['y', 'z'],      # S(y, z)
            ['z', 'u'],      # T(z, u)
            ['x', 'z', 'u'], # xz --> u
            ['y', 'u', 'x'], # yu --> x
        ];
        weights = [
            1.0,
            1.0,
            1.0,
            Inf,
            Inf,
        ]
    )
    @show(H)
    pb = polymatroid_bound(H)
    @assert pb ≈ 2.0
    fds = FD{Char}[
        FD(['x', 'z'], ['u']),
        FD(['y', 'u'], ['x']),
    ]
    pb = polymatroid_bound(H; fds = fds)
    @assert pb ≈ 1.5
end

# A variant of `test_polymatroid_bound1()`
function test_polymatroid_bound1_variant()
    println(repeat("=", 80))
    H = Hypergraph(
        ['x', 'y', 'z', 'u'],
        [
            ['x', 'y'],      # R(x, y)
            ['y', 'z'],      # S(y, z)
            ['z', 'u'],      # T(z, u)
            ['x', 'z', 'u'], # xz --> u
            ['y', 'u', 'x'], # yu --> x
        ];
        weights = [
            10.0,
            11.0,
            12.0,
            Inf,
            Inf,
        ]
    )
    @show(H)
    pb = polymatroid_bound(H)
    @assert pb ≈ 22.0 # 10 + 12
    fds = FD{Char}[
        FD(['x', 'z'], ['u']),
    ]
    dcs = DC{Char}[
        DC(['y', 'u'], ['x'], 2.0),
    ]
    pb = polymatroid_bound(H; fds = fds, dcs = dcs)
    @assert pb ≈ 35/2 # (10 + 11 + 12 + 2 + 0) / 2
end

# Q :- R(A, B), S(B, C), T(C, D), U(D, A)
# `U(D, A)` is an infinite relation but has two DCs from `A' to `D` and from `D` to `A`
function test_polymatroid_bound2()
    println(repeat("=", 80))
    H = Hypergraph(
        ['A', 'B', 'C', 'D'],
        [
            ['A', 'B'],      # R(A, B)
            ['B', 'C'],      # S(B, C)
            ['C', 'D'],      # T(C, D)
            ['D', 'A'],      # U(D, A)
        ];
        weights = [
            10.0,
            11.0,
            12.0,
            Inf,
        ]
    )
    @show(H)
    pb = polymatroid_bound(H)
    @assert pb ≈ 22.0 # 10 + 12
    fds = FD{Char}[
        FD(['D'], ['A']),
    ]
    dcs = DC{Char}[
        DC(['A'], ['D'], 3.0),
    ]
    pb = polymatroid_bound(H; fds = fds, dcs = dcs)
    @assert pb ≈ 36/2  # (10 + 11 + 12 + 3 + 0) / 2
end

# Run all tests
function test_all()
    test_4cycle()
    test_4cycle_with_fds()
    test_5cycle()
    test_5cycle_with_fds()
    test_6cycle()
    test_6cycle_with_fds()
    test_example_6()
    test_polymatroid_bound1()
    test_polymatroid_bound1_variant()
    test_polymatroid_bound2()
end

end
