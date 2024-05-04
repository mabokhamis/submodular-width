# import Pkg
# Pkg.add("JuMP")
# Pkg.add("Clp")
# Pkg.add("MathOptInterface")
# Pkg.add("Combinatorics")
# Pkg.add("DataStructures")

module OmegaSubmodularWidth

using JuMP
using Clp
using MathOptInterface
using Combinatorics
using DataStructures

"""
    Hypergraph{T}

Represents a hypergraph `H` whose vertices have type `T`. This struct has the following
fields:
    - `vars`: The set of vertices of `H`
    - `edges`: The set of hyperedges of `H`, each of which is a set of vertices
"""
mutable struct Hypergraph{T}
    vars::Vector{T}
    edges::Vector{Set{T}}

    # `var_index` maps a vertex `vars[i]` in `vars` to its index `i`
    var_index::Dict{T, Int}

    function Hypergraph(
        vars::Vector{T},
        edges::Vector{Vector{T}};
    ) where T
        @assert length(unique(vars)) == length(vars)
        var_index = Dict{T, Int}(var => i for (i, var) in enumerate(vars))
        @assert all(e ⊆ vars for e in edges)
        edges = map(edge -> Set{T}(edge), edges)
        return new{T}(vars, edges, var_index)
    end
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

function omega_submodular_width(H::Hypergraph{T}; verbose::Bool = true) where T
    @warn "This method does not handle the general case yet"

    n = length(H.vars)
    N = 2 ^ n

    function f(z::Int)::String where T
        return "h(" * join(map(string, sort(collect(unzip(H, z))))) * ")"
    end

    # initialize a linear program (LP)
    model = Model(Clp.Optimizer)
    set_optimizer_attribute(model, "LogLevel", 0)

    # Let `V` be the set of vertices of `H`. For each subset `U ⊆ V`, the LP contains a
    # corresponding variable `h[U]`
    @variable(model, h[0:N-1])

    # The LP contains the constraint `h[∅] = 0`
    verbose && println("\nZero Constraint:")
    @constraint(model, h[0] == 0.0)
    verbose && println("$(f(0)) = 0.0")

    # For each `X ⊆ Y ⊆ V`, the LP contains a constraint `h[X] ≤ h[Y]`. These are called
    # "monotonicity constraints"
    verbose && println("\n(Basic) Monotonicity Constraints:")
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
    verbose && println("\n(Basic) Submodularity Constraints:")
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
        E = zip(H, edge)
        @constraint(model, h[E] ≤ 1.0)
        verbose && println("$(f(E)) ≤ 1.0")
    end

    A = zip(H, Set([H.vars[1]]))
    B = zip(H, Set([H.vars[2]]))
    C = zip(H, Set([H.vars[3]]))
    ABC = zip(H, Set([H.vars[1], H.vars[2], H.vars[3]]))

    @variable(model, w >= 0.0)

    verbose && println("\nObjective Constraints:")
    @constraint(model, w <= h[ABC])
    verbose && println("w <= $(f(ABC))")
    @constraint(model, w <= h[A] + h[B])
    verbose && println("w <= $(f(A))+$(f(B))")

    # Finally, we set the objective of the LP to maximize `W`
    verbose && println("\nObjective:")
    @objective(model, Max, w)
    verbose && println("maximize w")

    optimize!(model)
    @assert termination_status(model) == MathOptInterface.OPTIMAL

    obj = objective_value(model)

    verbose && println("\nOptimal Solution:")
    sol = value.(h)
    for i = 0:N-1
        verbose && println("$(f(i)) = $(sol[i])")
    end

    return obj
end

H = Hypergraph(
    ["A", "B", "C"],
    [["A", "B"], ["B", "C"], ["A", "C"]]
)

println(omega_submodular_width(H))

end
