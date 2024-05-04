# import Pkg
# Pkg.add("JuMP")
# Pkg.add("Clp")
# Pkg.add("MathOptInterface")
# Pkg.add("Combinatorics")
# Pkg.add("DataStructures")
# Pkg.add("AutoHashEquals")

module OmegaSubmodularWidth

using JuMP
using Clp
using MathOptInterface
using Combinatorics
using DataStructures
using AutoHashEquals

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

"""
    name(U)

Given a set `U` of vertices, return a string representing `h(U)`
"""
function name(U::Set{T})::String where T
    return "h(" * join(map(string, sort(collect(U)))) * ")"
end

function omega_submodular_width(H::Hypergraph{T}; verbose::Bool = true) where T
    @warn "This method does not handle the general case yet"

    n = length(H.vars)
    N = 2 ^ n

    f(z::Int) = name(unzip(H, z))

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
    for edge in H.edges
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

# H = Hypergraph(
#     ["A", "B", "C"],
#     [["A", "B"], ["B", "C"], ["A", "C"]]
# )

# println(omega_submodular_width(H))

@auto_hash_equals struct Constant
    value::Float64
    symbol::String

    Constant(value::Float64, symbol::String = string(value)) = new(value, symbol)
end

Base.show(io::IO, c::Constant) = print(io, c.symbol)

abstract type Term{T} end

@auto_hash_equals struct Sum{T} <: Term{T}
    args::Dict{Set{T},Constant}
end

Sum(args::Dict{Set{T},Constant}) where T = Sum{T}(args)

@auto_hash_equals struct Min{T} <: Term{T}
    args::Vector{Term{T}}
end

Min(args::Vector{<:Term{T}}) where T = Min{T}(args)

@auto_hash_equals struct Max{T} <: Term{T}
    args::Vector{Term{T}}
end

Max(args::Vector{<:Term{T}}) where T = Max{T}(args)

function pretty_print(io::IO, s::Sum; indent::Int = 0)
    margin = repeat(" ", indent)
    t = SortedDict(name(x) => c for (x, c) in s.args)
    print(io, margin, join(("$c*$x" for (x, c) in t), " + "))
end

function pretty_print(io::IO, m::Union{Min, Max}; indent::Int = 0)
    margin = repeat(" ", indent)
    t = m isa Min ? "Min" : "Max"
    println(io, "$(margin)$t{")
    for arg ∈ m.args
        pretty_print(io, arg; indent = indent + 4)
        println(io, ",")
    end
    print(io, "$(margin)}")
end

function Base.show(io::IO, t::Term)
    pretty_print(io, t)
end

_same_type(x::T, y::T) where T = true
_same_type(x, y) = false

function coefficient(sum::Sum{T}, x::Set{T}) where T
    return  haskey(sum.args, x) ? sum.args[x].value : 0.0
end

function Base.:(<=)(a::Sum, b::Sum)
    return all(coefficient(a, x) <= coefficient(b, x) for x ∈ keys(a.args) ∪ keys(b.args))
end

function Base.:(<=)(a::Sum, b::Min)
    return all(a ≤ arg for arg ∈ b.args)
end

function Base.:(<=)(a::Sum, b::Max)
    return any(a ≤ arg for arg ∈ b.args)
end

function Base.:(<=)(a::Min, b::Term)
    return any(arg ≤ b for arg ∈ a.args)
end

function Base.:(<=)(a::Max, b::Term)
    return all(arg ≤ b for arg ∈ a.args)
end

_min_subsumed_by(i, a, j, b) = b <= a && (!(a <= b) || a <= b && j < i)
_max_subsumed_by(i, a, j, b) = a <= b && (!(b <= a) || b <= a && j < i)

function minimal_args(args::Vector{<:Term{T}}; subsumed_by::Function = _min_subsumed_by) where {T}
    new_args = Vector{Term{T}}()
    for (i, a) ∈ enumerate(args)
        if any(subsumed_by(i, a, j, b) for (j, b) ∈ enumerate(args) if j != i)
            continue
        end
        push!(new_args, a)
    end
    return new_args
end

_simplify(s::Sum) = (false, s)

function _simplify(m::Union{Min{T},Max{T}}) where T
    new_args = unique(m.args)
    if length(new_args) != length(m.args)
        return (true, typeof(m)(new_args))
    end

    if any(_same_type(arg, m) for arg ∈ m.args)
        new_args = Vector{Term{T}}()
        for arg ∈ m.args
            if _same_type(arg, m)
                append!(new_args, arg.args)
            else
                push!(new_args, arg)
            end
        end
        return (true, typeof(m)(new_args))
    end

    if m isa Min || m isa Max
        new_args = minimal_args(m.args;
            subsumed_by = m isa Min ? _min_subsumed_by : _max_subsumed_by)
        if length(new_args) != length(m.args)
            return (true, typeof(m)(new_args))
        end
    end

    if m isa Min && any(arg isa Max for arg ∈ m.args)
        new_args = Vector{Vector{Term{T}}}()
        for arg ∈ m.args
            if arg isa Max
                push!(new_args, arg.args)
            else
                push!(new_args, [arg])
            end
        end
        selectors = Iterators.product(new_args...,)
        new_args = Vector{Term{T}}()
        for β ∈ selectors
            push!(new_args, Min(collect(β)))
        end
        return (true, Max(new_args))
    end

    return (false, m)
end

function simplify(t::Term{T}) where T
    if t isa Min || t isa Max
        t = typeof(t)([simplify(arg) for arg ∈ t.args])
    end
    b = true
    while b
        (b, t) = _simplify(t)
    end
    if t isa Min || t isa Max
        t = typeof(t)([simplify(arg) for arg ∈ t.args])
    end
    return t
end

exp =
Max([
    Min([
        Max([
            Sum(Dict(Set(["A"]) => Constant(1.0))),
            Sum(Dict(Set(["B"]) => Constant(1.0))),
        ]),
        Max([
            Sum(Dict(Set(["C"]) => Constant(1.0))),
            Sum(Dict(Set(["D"]) => Constant(1.0))),
        ]),
        Max([
            Sum(Dict(Set(["C"]) => Constant(1.0))),
            Sum(Dict(Set(["D"]) => Constant(1.0))),
        ]),
    ]),
    Sum(Dict(Set(["A"]) => Constant(1.0), Set(["B"]) => Constant(0.9))),
    Sum(Dict(Set(["A"]) => Constant(0.9), Set(["B"]) => Constant(1.9))),
])

println(exp)

exp = simplify(exp)

println(exp)

end
