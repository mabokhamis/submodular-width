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
        edges::Vector{<:Union{Vector{T},Set{T}}};
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
    @assert U ⊆ H.vars
    return zip(H.var_index, U)
end

function zip(var_index::Dict{T, Int}, U::Set{T})::Int where T
    z = 0
    for x in U
        z |= (1 << (var_index[x] - 1))
    end
    return z
end

"""
    unzip(H, z)

Given a hypergraph `H` and an integer `z` representing a subset `U` of the vertices of `H`
(that was encoded using `z = zip(H, U)`), return `U`
"""
function unzip(H::Hypergraph{T}, z::Int)::Set{T} where T
    return unzip(H.vars, z)
end

function unzip(vars::Vector{T}, z::Int)::Set{T} where T
    set = Set{T}()
    i = 1
    while z != 0
        if z & 1 == 1
            push!(set, vars[i])
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

struct Constant
    value::Float64
    symbol::String

    Constant(value::Float64, symbol::String = string(value)) = new(value, symbol)
end

Base.:(==)(x::Constant, y::Constant) = x.value == y.value
Base.hash(x::Constant, h::UInt) = hash(x.value, h)

Base.show(io::IO, c::Constant) = print(io, c.symbol)

abstract type Term{T} end

@auto_hash_equals struct Sum{T} <: Term{T}
    args::Dict{Set{T},Constant}
end

function Sum(args::Dict{Set{T},Constant}) where T
    args = Dict(v => c for (v, c) in args if !isempty(v) && c.value != 0.0)
    return Sum{T}(args)
end

function coefficient(sum::Sum{T}, x::Set{T}) where T
    return  haskey(sum.args, x) ? sum.args[x].value : 0.0
end

function Base.:(-)(a::Sum{T}, b::Sum{T}) where T
    args = Dict{Set{T}, Constant}()
    for x ∈ keys(a.args) ∪ keys(b.args)
        c = coefficient(a, x) - coefficient(b, x)
        c != 0.0 && (args[x] = Constant(c))
    end
    return Sum(args)
end

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

function Base.:(<=)(a::Sum, b::Sum)
    return is_non_negative(b - a)
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
        new_args = [simplify(arg) for arg ∈ t.args]
        t = typeof(t)(new_args)
    end
    rewritten = false
    b = true
    while b
        (b, t) = _simplify(t)
        rewritten |= b
    end
    return rewritten ? simplify(t) : t
end

# exp =
# Max([
#     Min([
#         Max([
#             Sum(Dict(Set(["A"]) => Constant(1.0))),
#             Sum(Dict(Set(["B"]) => Constant(1.0))),
#         ]),
#         Max([
#             Sum(Dict(Set(["C"]) => Constant(1.0))),
#             Sum(Dict(Set(["D"]) => Constant(1.0))),
#         ]),
#         Max([
#             Sum(Dict(Set(["C"]) => Constant(1.0))),
#             Sum(Dict(Set(["D"]) => Constant(1.0))),
#         ]),
#     ]),
#     Sum(Dict(Set(["A"]) => Constant(1.0), Set(["B"]) => Constant(0.9))),
#     Sum(Dict(Set(["A"]) => Constant(0.9), Set(["B"]) => Constant(1.9))),
# ])

# println(exp)

# @time exp = simplify(exp)

# println(exp)

function create_matrix_multiplication(
    X::Set{T},
    Y::Set{T},
    Z::Set{T},
    W::Set{T},
    ω::Float64
) where T
    @assert isempty(X ∩ Y) && isempty(X ∩ Z) && isempty(X ∩ W) && isempty(Y ∩ Z) &&
        isempty(Y ∩ W) && isempty(Z ∩ W)
    one = Constant(1.0)
    ω2 = Constant(ω-2, "ω'")
    return Max([
        Sum(Dict(
            X => one,
            Y => one,
            Z => ω2,
            W => one,
        )),
        Sum(Dict(
            X => one,
            Y => ω2,
            Z => one,
            W => one,
        )),
        Sum(Dict(
            X => ω2,
            Y => one,
            Z => one,
            W => one,
        )),
    ])
end

# e = create_matrix_multiplication(Set(["A"]), Set(["B"]), Set(["C"]), Set{String}(), 2.5)
# println(simplify(e))

function eliminate_variable(H::Hypergraph{T}, x::T, ω::Float64) where T
    Nx = [E for E ∈ H.edges if x ∈ E]
    Px = [E for E ∈ H.edges if x ∉ E]
    U = reduce(union!, Nx; init = Set{T}())
    P = copy(U)
    for E ∈ Px
        setdiff!(P, E)
    end
    args = Vector{Term{T}}()
    push!(args, Sum(Dict(U => Constant(1.0))))
    k = length(Nx)
    for s ∈ powerset(collect(1:k))
        (isempty(s) || length(s) == k) && continue
        A = reduce(union!, Nx[i] for i ∈ s; init = Set{T}())
        B = reduce(union!, Nx[i] for i ∈ 1:k if i ∉ s; init = Set{T}())
        Y = A ∩ B ∩ P
        W = setdiff(A ∩ B, P)
        X = setdiff(A, B)
        Z = setdiff(B, A)
        (isempty(X) || isempty(Y) || isempty(Z)) && continue
        arg = create_matrix_multiplication(X, Y, Z, W, ω)
        push!(args, arg)
    end
    new_vars = setdiff(H.vars, P)
    new_edges = [[setdiff(E, P) for E ∈ Px]; setdiff(U, P)]
    filter(E -> !isempty(E), new_edges)
    new_H = Hypergraph(new_vars, new_edges)
    expr = simplify(Min(args))
    return new_H, expr
end

function eliminate_variables(H::Hypergraph{T}, ω::Float64) where T
    min_args = Vector{Term{T}}()
    for π ∈ permutations(H.vars)
        new_H = H
        max_args = Vector{Term{T}}()
        for x ∈ π
            x ∈ new_H.vars || continue
            (new_H, expr) = eliminate_variable(new_H, x, ω)
            push!(max_args, expr)
        end
        arg = Max(max_args)
        push!(min_args, arg)
    end
    return simplify(Min(min_args))
end

function solve(A, b)
    model = Model(Clp.Optimizer)
    set_optimizer_attribute(model, "LogLevel", 0)
    @variable(model, x[1:size(A, 2)] >= 0)
    b1 = b .- 1e-7
    b2 = b .+ 1e-7
    @constraint(model, A * x >= b1)
    @constraint(model, A * x <= b2)
    optimize!(model)
    return termination_status(model) == MathOptInterface.OPTIMAL
end

function is_non_negative(sum::Sum{T}) where T
    vars = collect(reduce(union!, keys(sum.args); init = Set{T}()))
    var_index = Dict{T, Int}(var => i for (i, var) in enumerate(vars))
    n = length(vars)
    N = 2 ^ n

    A = Vector{Vector{Float64}}()
    c = zeros(N)
    c[0 + 1] = 1.0
    push!(A, c)
    c = zeros(N)
    c[0 + 1] = -1.0
    push!(A, c)

    for y = 0:n-1
        Y = N - 1
        X = Y & ~(1 << y)
        c = zeros(N)
        c[Y + 1] = 1.0
        c[X + 1] = -1.0
        push!(A, c)
    end

    for X = 0:N-1, y = 0:n-1, z = y+1:n-1
        if (X & (1 << y) == 0) && (X & (1 << z) == 0)
            Y = X | (1 << y)
            Z = X | (1 << z)
            W = Y | (1 << z)
            c = zeros(N)
            c[Y + 1] = 1.0
            c[Z + 1] = 1.0
            c[X + 1] = -1.0
            c[W + 1] = -1.0
            push!(A, c)
        end
    end

    A = hcat(A...)

    b = zeros(N)
    for (v, c) ∈ sum.args
        b[zip(var_index, v) + 1] = c.value
    end

    return solve(A, b)
end

# s = Sum(Dict(
#     Set(["A", "B", "C"]) => Constant(-1.0),
#     Set(["A", "B"]) => Constant(0.5),
#     Set(["B", "C"]) => Constant(0.5),
#     Set(["C", "A"]) => Constant(0.5),
# ))

# println(is_non_negative(s))

function omega_submodular_width(H::Hypergraph{T}, m::Min{T}; verbose::Bool = true) where T

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

    @variable(model, w >= 0.0)

    verbose && println("\nObjective Constraints:")
    for s ∈ m.args
        @constraint(model, w ≤ sum(c.value * h[zip(H, X)] for (X, c) ∈ s.args))
        verbose && println("w ≤ $(join(["$(c.value) * $(f(zip(H, X)))" for (X, c) ∈ s.args], " + "))")
    end

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

# m = Min([
#     Sum(Dict(Set(["A", "B", "C"]) => Constant(1.0))),
#     Sum(Dict(Set(["A"]) => Constant(1.0), Set(["B"]) => Constant(1.0))),
# ])

# println(omega_submodular_width(H, m; verbose = true))

function omega_submodular_width(H::Hypergraph{T}, ω::Float64; verbose::Bool = true) where T
    expr = eliminate_variables(H, ω)
    println(expr)
    @assert expr isa Max
    width = 0.0
    for arg ∈ expr.args
        @assert arg isa Min
        width = max(width, omega_submodular_width(H, arg; verbose = verbose))
    end
    return width
end

H = Hypergraph(
    ["A", "B", "C"],
    [["A", "B"], ["B", "C"], ["A", "C"]]
)

ω = 2.5
w = omega_submodular_width(H, ω; verbose = false)
println(w)
println(2 * ω / (ω + 1))

#-----------------------------------------------

# H = Hypergraph(
#     ["A", "B1", "B2", "B3"],
#     [["B1", "B2", "B3"], ["A", "B1"], ["A", "B2"], ["A", "B3"]]
# )

# ω = 3.0
# w = omega_submodular_width(H, ω; verbose = false)
# println(w)
# println(1 + 2 * ω / (2ω + 3))

# expr = Max([
#     Min([
#         Sum(Dict(Set(["A", "B", "C"]) => Constant(1.0))),
#         Sum(Dict(Set(["A"]) => Constant(1.0), Set(["B"]) => Constant(1.0))),
#         Sum(Dict(Set(["B"]) => Constant(1.0), Set(["C"]) => Constant(1.0))),
#     ]),
#     Min([
#         Sum(Dict(Set(["A", "B", "C"]) => Constant(1.0))),
#         Sum(Dict(Set(["B"]) => Constant(1.0), Set(["C"]) => Constant(1.0))),
#         Sum(Dict(Set(["A"]) => Constant(1.0), Set(["B"]) => Constant(1.0))),
#     ]),
#     Min([
#         Sum(Dict(Set(["A", "B", "C"]) => Constant(1.0))),
#         Sum(Dict(Set(["B"]) => Constant(1.0), Set(["C"]) => Constant(1.0))),
#         Sum(Dict(Set(["A"]) => Constant(1.0), Set(["B"]) => Constant(1.0))),
#         Sum(Dict(Set(["A"]) => Constant(1.0), Set(["C"]) => Constant(1.0))),
#     ]),
# ])

# expr2 = simplify(expr)

# println(repeat("%", 500))
# # println(expr)
# # println(repeat("-", 40))
# # println(expr2)
# x = expr2.args[1]
# y = expr2.args[3]
# println(x)
# println(y)
# println(x <= y)

end
