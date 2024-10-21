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

#-------------------------------------------------------------------------------------------

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
end

function Hypergraph(
    vars::Vector{T},
    edges::Vector{<:Union{Vector{T},Set{T}}}
) where T
    @assert length(unique(vars)) == length(vars)
    var_index = Dict{T, Int}(var => i for (i, var) in enumerate(vars))
    @assert all(e ⊆ vars for e in edges)
    edges = map(edge -> Set{T}(edge), edges)
    @assert all(reduce(union!, edges; init = Set{T}()) == Set{T}(vars))
    return Hypergraph{T}(vars, edges, var_index)
end

# Copy a hypergraph
function Base.copy(H::Hypergraph{T}) where T
    return Hypergraph{T}(
        copy(H.vars),
        [copy(E) for E in H.edges],
        copy(H.var_index)
    )
end

function ∂(H::Hypergraph{T}, X::Set{T})::Vector{Set{T}} where T
    return filter(E -> !isempty(E ∩ X), H.edges)
end

function U(H::Hypergraph{T}, X::Set{T})::Set{T} where T
    return reduce(union!, ∂(H, X); init = Set{T}())
end

function N(H::Hypergraph{T}, X::Set{T})::Set{T} where T
    return setdiff(U(H, X), X)
end

function eliminate!(H::Hypergraph{T}, X::Set{T})::Hypergraph{T} where T
    new_E = N(H, X)
    filter!(E -> isempty(E ∩ X), H.edges)
    push!(H.edges, new_E)
    setdiff!(H.vars, X)
    return H
end

"""
    zip(H, U)

Given a hypergraph `H` and a subset `U` of the vertices of `H`, encode `U` as a string of
bits. For example, `{v1, v3, v4, v8}` is encoded as the binary string `10001101`
"""
function zip(H::Hypergraph{T}, U::Union{Set{T},Vector{T}})::Int where T
    return zip(H.var_index, U)
end

function zip(var_index::Dict{T, Int}, U::Union{Set{T},Vector{T}})::Int where T
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
function name(U::Union{Set{T},Vector{T}})::String where T
    return "h(" * join(map(string, sort!(collect(U)))) * ")"
end

#-------------------------------------------------------------------------------------------

struct Constant
    value::Float64
    symbol::String

    Constant(value::Number, symbol::String = string(value)) = new(Float64(value), symbol)
end

Base.:(==)(x::Constant, y::Constant) = x.value == y.value
Base.hash(x::Constant, h::UInt) = hash(x.value, h)

Base.show(io::IO, c::Constant) = print(io, c.symbol)

abstract type Term{T} end

@auto_hash_equals struct Sum{T} <: Term{T}
    args::Dict{Set{T},Constant}

    function Sum(args::Dict{Set{T},Constant}, τ::Number = 1e-7) where T
        args = Dict(v => c for (v, c) in args if !isempty(v) && abs(c.value) > τ)
        return new{T}(args)
    end
end

@auto_hash_equals struct Min{T} <: Term{T}
    args::Set{Term{T}}

    function Min(args::Union{Set{<:Term{T}},Vector{<:Term{T}}}) where T
        new_args = Set{Term{T}}()
        for arg ∈ args
            if arg isa Min || arg isa Max && length(arg.args) == 1
                union!(new_args, arg.args)
            else
                push!(new_args, arg)
            end
        end
        return new{T}(new_args)
    end
end

@auto_hash_equals struct Max{T} <: Term{T}
    args::Set{Term{T}}

    function Max(args::Union{Set{<:Term{T}},Vector{<:Term{T}}}) where T
        new_args = Set{Term{T}}()
        for arg ∈ args
            if arg isa Max || arg isa Min && length(arg.args) == 1
                union!(new_args, arg.args)
            else
                push!(new_args, arg)
            end
        end
        return new{T}(new_args)
    end
end

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

function _distribute(m::Union{Min{T},Max{T}}) where T
    if m isa Min && any(arg isa Max for arg ∈ m.args)
        new_args = Vector{Vector{Term{T}}}()
        num_args = 1
        args = collect(m.args)
        for arg ∈ args
            if arg isa Max
                push!(new_args, collect(arg.args))
                num_args *= length(arg.args)
            else
                push!(new_args, [arg])
            end
        end
        num_args >= 15 && println("Begin _distribute $num_args")
        selectors = Iterators.product(new_args...,)
        new_args = Vector{Term{T}}()
        for β ∈ selectors
            push!(new_args, Min(collect(β)))
        end
        num_args >= 15 && println("End _distribute $num_args")
        return Max(new_args)
    end
    return m
end

function distribute(m::Term)
    return bottomup(m, _distribute)
end

function bottomup(m::Term{T}, rewrite::Function) where T
    if m isa Min || m isa Max
        new_args = [bottomup(arg, rewrite) for arg ∈ m.args]
        m = m isa Min ? Min(new_args) : Max(new_args)
        return rewrite(m)
    end
    return m
end

function topdown(m::Term{T}, rewrite::Function) where T
    if m isa Min || m isa Max
        m = rewrite(m)
        new_args = [topdown(arg, rewrite) for arg ∈ m.args]
        m = m isa Min ? Min(new_args) : Max(new_args)
    end
    return m
end

function _distribute_and_simplify(m::Union{Min{T},Max{T}}) where T
    m = _distribute(m)
    num_args = length(m.args)
    if num_args <= 50
        num_args >= 10 && println("Begin _remove_subsumed_args $num_args")
        m = _remove_subsumed_args(m)
        num_args >= 10 && println("End _remove_subsumed_args $num_args")
    end
    return m
end

function distribute_and_simplify(m::Term)
    return bottomup(m, _distribute_and_simplify)
end

#-------------------------------------------------------------------------------------------

function coefficient(sum::Sum{T}, x::Set{T}) where T
    return  haskey(sum.args, x) ? sum.args[x].value : 0.0
end

function Base.:(-)(a::Sum{T}, b::Sum{T}) where T
    args = Dict{Set{T}, Constant}()
    for x ∈ keys(a.args) ∪ keys(b.args)
        args[x] = Constant(coefficient(a, x) - coefficient(b, x))
    end
    return Sum(args)
end

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

function Base.:(<=)(a::Min, b::Min)
    return all(any(a2 <= b2 for a2 ∈ a.args) for b2 ∈ b.args)
end

function Base.:(<=)(a::Max, b::Max)
    return all(any(a2 <= b2 for b2 ∈ b.args) for a2 ∈ a.args)
end

_min_subsumed_by(i, a, j, b) = b <= a && (!(a <= b) || j < i)
_max_subsumed_by(i, a, j, b) = a <= b && (!(b <= a) || j < i)

function _minimal_args(args::Vector{<:Term{T}}; subsumed_by::Function = _min_subsumed_by) where {T}
    new_args = Vector{Term{T}}()
    for (i, a) ∈ enumerate(args)
        any(subsumed_by(i, a, j, b) for (j, b) ∈ enumerate(args) if j != i) && continue
        push!(new_args, a)
    end
    return new_args
end

function _remove_subsumed_args(m::Min{T}) where T
    new_args = _minimal_args(collect(m.args); subsumed_by = _min_subsumed_by)
    return Min(new_args)
end

function _remove_subsumed_args(m::Max{T}) where T
    new_args = _minimal_args(collect(m.args); subsumed_by = _max_subsumed_by)
    return Max(new_args)
end

function remove_subsumed_args(m::Term)
    return bottomup(m, _remove_subsumed_args)
end

#-------------------------------------------------------------------------------------------

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

#-------------------------------------------------------------------------------------------

function MM(
    X::Set{T},
    Y::Set{T},
    Z::Set{T},
    G::Set{T},
    ω::Number,
    verbose::Bool = false
) where T
    @assert isempty(X ∩ Y) && isempty(X ∩ Z) && isempty(X ∩ G) && isempty(Y ∩ Z) &&
        isempty(Y ∩ G) && isempty(Z ∩ G)
    verbose && println("MM($X ; $Y ; $Z | $G)")
    one = Constant(1.0, "1")
    γ = Constant(ω-2, "γ")
    α = Constant(-ω+1, "α")
    return Max([
        Sum(Dict(
            X ∪ G => one,
            Y ∪ G => one,
            Z ∪ G => γ,
            G => α,
        )),
        Sum(Dict(
            X ∪ G => one,
            Y ∪ G => γ,
            Z ∪ G => one,
            G => α,
        )),
        Sum(Dict(
            X ∪ G => γ,
            Y ∪ G => one,
            Z ∪ G => one,
            G => α,
        )),
    ])
end

function MM(
    X::Vector{T}, Y::Vector{T}, Z::Vector{T}, G::Vector{T}, ω::Number, verbose::Bool = false
) where T
    return MM(Set{T}(X), Set{T}(Y), Set{T}(Z), Set{T}(G), ω, verbose)
end

function min_U_EMM(H::Hypergraph{T}, X::Set{T}, ω::Number, verbose::Bool = false) where T
    isempty(X) && return nothing

    ∂_X = ∂(H, X)
    U_X = U(H, X)
    any(E == U_X for E ∈ ∂_X) && return nothing

    args = Vector{Term{T}}()
    push!(args, Sum(Dict(U_X => Constant(1.0, "1"))))
    comb = Set{Tuple{Set{T},Set{T},Set{T}}}()
    k = length(∂_X)
    for AA ∈ powerset(collect(1:k))
        (isempty(AA) || length(AA) == k) && continue
        for AB ∈ powerset(AA)
            length(AB) == length(AA) && continue
            A = reduce(union!, ∂_X[i] for i ∈ AA; init = Set{T}())
            B = reduce(union!, ∂_X[i] for i ∈ 1:k if i ∉ AA || i ∈ AB; init = Set{T}())
            A_cap_B = A ∩ B
            X ⊆ A_cap_B || continue
            A_diff_B = setdiff(A, B)
            isempty(A_diff_B) && continue
            B_diff_A = setdiff(B, A)
            isempty(B_diff_A) && continue
            G1 = setdiff(A_cap_B, X)
            for G2 ∈ powerset(collect(setdiff(A ∪ B, A_cap_B)))
                G = G1 ∪ G2
                Y = setdiff(A_diff_B, G)
                isempty(Y) && continue
                Z = setdiff(B_diff_A, G)
                isempty(Z) && continue
                (Y, Z, G) ∈ comb && continue
                (Z, Y, G) ∈ comb && continue
                push!(comb, (Y, Z, G))
                arg = MM(X, Y, Z, G, ω, verbose)
                push!(args, arg)
            end
        end
    end
    return Min(args)
end

function min_elimination_cost(H::Hypergraph{T}, ω::Number, verbose::Bool = false) where T
    isempty(H.vars) && return nothing

    args = Vector{Term{T}}()
    for X ∈ powerset(H.vars)
        isempty(X) && continue
        expr1 = min_U_EMM(H, Set(X), ω, verbose)
        H2 = copy(H)
        eliminate!(H2, Set(X))
        expr2 = min_elimination_cost(H2, ω, verbose)
        arg = if isnothing(expr1) || isnothing(expr2)
            isnothing(expr1) ? expr2 : expr1
        else
            Max([expr1, expr2])
        end
        isnothing(arg) && return nothing
        push!(args, arg)
    end
    return Min(args)
end

#-------------------------------------------------------------------------------------------

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
    Y = N - 1
    for y = 0:n-1
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
        @assert s isa Sum
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

    polymatroid = Dict{Set{T}, Float64}()
    verbose && println("\nOptimal Primal Solution:")
    sol = value.(h)
    for i = 0:N-1
        polymatroid[unzip(H, i)] = sol[i]
        verbose && println("$(f(i)) = $(sol[i])")
    end

    return (obj, polymatroid)
end

function omega_submodular_width(H::Hypergraph{T}, ω::Number; verbose::Bool = true, expr = nothing) where T
    isnothing(expr) && (expr = min_elimination_cost(H, ω, verbose))
    println(expr)
    expr = distribute(expr)
    if !(expr isa Max)
        expr = Max([expr])
    end
    width = 0.0
    witness = nothing
    for (i, arg) ∈ enumerate(expr.args)
        (i%10 == 1) && println("C: $i of $(length(expr.args))")
        if !(arg isa Min)
            arg = Min([arg])
        end
        bound, h = omega_submodular_width(H, arg; verbose = false)
        bound2 = eval(arg, h)
        @assert abs(bound - bound2) < 1e-6 """
         - bound = $bound
         - bound2 = $bound2
        """
        if bound > width
            width = bound
            witness = h
        end
    end
    return (width, witness)
end

#-------------------------------------------------------------------------------------------

function eval(c::Constant, h::Dict{Set{T}, Float64}) where T
    return c.value
end

function eval(s::Sum{T}, h::Dict{Set{T}, Float64}) where T
    return sum(c.value * h[x] for (x, c) ∈ s.args; init = 0.0)
end

function eval(m::Min{T}, h::Dict{Set{T}, Float64}) where T
    return minimum(eval(arg, h) for arg ∈ m.args; init = Inf)
end

function eval(m::Max{T}, h::Dict{Set{T}, Float64}) where T
    return maximum(eval(arg, h) for arg ∈ m.args; init = -Inf)
end

function Base.show(io::IO, h::Dict{Set{T}, Float64}) where T
    v = [(x, v) for (x, v) ∈ h]
    sort!(v, by = x -> (length(x[1]), name(x[1])))
    for (x, v) ∈ v
        println(io, "$(name(x)) = $v")
    end
end

#-------------------------------------------------------------------------------------------

function test1()
    ω = 2.5
    t = Min([
        MM(["X"], ["Y"], ["Z"], ["W"], ω),
        MM(["Y"], ["Z"], ["X"], ["W"], ω),
        MM(["Z"], ["X"], ["Y"], ["W"], ω),
        MM(["Z", "T"], ["X"], ["Y"], ["W"], ω),
        MM(["T"], ["X"], ["Y"], ["W"], ω),
    ])
    println(t)
    t = remove_subsumed_args(t)
    println(t)
    t = distribute(t)
    println(t)
end

# test1()

function test2()
    ω = 2.5
    H = Hypergraph(
        ["A", "B", "C"],
        [["A", "B"], ["B", "C"], ["A", "C"]]
    )
    # X = ["A"]
    # t = min_U_EMM(H, Set(X), ω, true)
    # println(t)
    w = omega_submodular_width(H, ω; verbose = true)
    println(w)
    println(2 * ω / (ω + 1))
end

# test2()

function test3()
    ω = 2.5
    H = Hypergraph(
        ["A", "B", "C", "D"],
        [["A", "B"], ["A", "C"], ["A", "D"], ["B", "C"], ["B", "D"], ["C", "D"]]
    )
    # X = ["A"]
    # t = min_U_EMM(H, Set(X), ω, true)
    # println(t)
    w = omega_submodular_width(H, ω; verbose = true)
    println(w)
    println((ω + 1)/2)
end

# test3()

#-------------------------------------------------------------------------------------------

# H = Hypergraph(
#     ["A", "B", "C"],
#     [["A", "B"], ["B", "C"], ["A", "C"]]
# )

# ω = 2
# w = omega_submodular_width(H, ω; verbose = false)
# println(w)
# println(2 * ω / (ω + 1))

#-----------------------------------------------

# H = Hypergraph(
#     ["A", "B1", "B2", "B3"],
#     [["B1", "B2", "B3"], ["A", "B1"], ["A", "B2"], ["A", "B3"]]
# )

# ω = 2
# w = omega_submodular_width(H, ω; verbose = false)
# println(w)
# println(1 + 2 * ω / (2ω + 3))

#=
for ω = 3.0
    omega_submodular_width = 1.666666666666667
    1 + 2 * ω / (2ω + 3)   = 1.6666666666666665
for ω = 2.5
    omega_submodular_width = 1.6
    1 + 2 * ω / (2ω + 3)   = 1.625
for ω = 2.0
    omega_submodular_width = 1.5
    1 + 2 * ω / (2ω + 3)   = 1.5714285714285714
=#

#-----------------------------------------------

# H = Hypergraph(
#     ["A", "B", "C", "D"],
#     [["A", "B"], ["A", "C"], ["A", "D"], ["B", "C"], ["B", "D"], ["C", "D"]]
# )

# ω = 2.0
# w = omega_submodular_width(H, ω; verbose = false)
# println(w)
# println((ω + 1)/2)

#=
for ω = 3.0
    omega_submodular_width = 2
    (ω + 1)/2   = 2
for ω = 2.5
    omega_submodular_width = 1.8571428571428577
    (ω + 1)/2   = 1.75
for ω = 2.0
    omega_submodular_width = 1.6666666666666
    (ω + 1)/2   = 1.5
=#

#-----------------------------------------------

# # 4-cycle:
# H = Hypergraph(
#     ["A", "B", "C", "D"],
#     [["A", "B"], ["B", "C"], ["C", "D"], ["D", "A"]]
# )

# ω = 2
# w = omega_submodular_width(H, ω; verbose = false)
# println(w)
# #=
# for ω = 3.0
#     omega_submodular_width = 1.5
# for ω = 2.0
#     omega_submodular_width = 1.4
# =#

#-----------------------------------------------

# # 5-cycle:
# H = Hypergraph(
#     ["A", "B", "C", "D", "E"],
#     [["A", "B"], ["B", "C"], ["C", "D"], ["D", "E"], ["E", "A"]]
# )

# ω = 2
# w = omega_submodular_width(H, ω; verbose = false)
# println(w)
# #=
# for ω = 3.0
#     omega_submodular_width = 5/3
# for ω = 2.0
#     omega_submodular_width = 1.5
# =#

#-----------------------------------------------

# H = Hypergraph(
#     ["A", "B", "C", "D", "E"],
#     [["A", "B", "D"], ["A", "B", "E"], ["A", "C"], ["B", "C"], ["C", "D", "E"]]
# )

# ω = 2
# w = omega_submodular_width(H, ω; verbose = false)
# println(w)
# println(2 - 5/(6 * ω + 7))

#-----------------------------------------------

# H = Hypergraph(
#     ["A", "B", "C", "D", "E"],
#     [["A", "B", "C", "D"], ["B", "C", "E"], ["B", "D", "E"], ["C", "D", "E"], ["A", "E"]];
# )

# ω = 2
# w = omega_submodular_width(H, ω; verbose = false)
# println(w)                      # 1.4285714285714286
# println(2 - 10/(6 * ω + 7))     # 1.473684210526316

#-----------------------------------------------

# H = Hypergraph(
#     ["A", "B", "C", "D"],
#     [["A", "B", "C"], ["B", "C", "D"], ["C", "D", "A"], ["D", "A", "B"]]
# )

# ω = 2
# w = omega_submodular_width(H, ω; verbose = false)
# println(w)
# #=
# for ω = 3.0
#     omega_submodular_width = 5/3
# for ω = 2.0
#     omega_submodular_width = 1.5
# =#

#-----------------------------------------------

# H = Hypergraph(
#     ["A1", "A2", "B1", "B2", "C1", "C2"],
#     [["A1", "A2", "B1", "B2"], ["B1", "B2", "C1", "C2"], ["A1", "A2", "C1", "C2"]];
# )

# ω = 2
# w = omega_submodular_width(H, ω; verbose = false)
# println(w)

#-----------------------------------------------

# H = Hypergraph(
#     ["Y", "X1", "X2", "X3"],
#     [["X1", "Y"], ["X2", "Y"], ["X3", "Y"], ["X1", "X2", "X3"]]
# )

# ω = 2
# w = omega_submodular_width(H, ω; verbose = false)
# println(w)
# println(2 - 1/ω)

#-----------------------------------------------

# H = Hypergraph(
#     ["A", "B", "C", "D"],
#     [["A", "B"], ["B", "C"], ["C", "D"], ["D", "A"]]
# )

# ω = 2
# w = omega_submodular_width(H, ω; verbose = false)
# println(w)

#-----------------------------------------------

# H = Hypergraph(
#     ["A", "B", "C", "D", "E"],
#     [["A", "B"], ["B", "C"], ["C", "D"], ["D", "E"], ["E", "A"]]
# )

# ω = 2
# w = omega_submodular_width(H, ω; verbose = false)
# println(w)

#-----------------------------------------------

H = Hypergraph(
    ["Y", "X1", "X2", "X3", "X4"],
    [["X1", "Y"], ["X2", "Y"], ["X3", "Y"], ["X4", "Y"], ["X1", "X2", "X3", "X4"]]
)

ω = 2.5

full_expr = min_elimination_cost(H, ω, false)
println(full_expr)
println()
println()
println()

expr = Min([
    Sum(Dict(Set(["Y", "X1", "X2", "X3", "X4"]) => Constant(1.0))),
    MM(["X1", "X2"], ["X3", "X4"], ["Y"], String[], ω),
    MM(["X1", "X3"], ["X2", "X4"], ["Y"], String[], ω),
    MM(["X1", "X4"], ["X2", "X3"], ["Y"], String[], ω),
    MM(["X1"], ["X2", "X3", "X4"], ["Y"], String[], ω),
    MM(["X2"], ["X1", "X3", "X4"], ["Y"], String[], ω),
    MM(["X3"], ["X1", "X2", "X4"], ["Y"], String[], ω),
    MM(["X4"], ["X1", "X2", "X3"], ["Y"], String[], ω),
    MM(["X4"], ["X2", "X3"], ["X1", "Y"], String[], ω),
])
(w, h) = omega_submodular_width(H, ω; expr, verbose = false)
println(w)
full_w = eval(full_expr, h)
@assert abs(w - full_w) < 1e-6 """
    - w = $w
    - full_w = $full_w
"""
println(h)

#-------------------------------------------------------------------------------------------

function lower_bound_triangle_tree(k::Int, ω::Number, ϵ::Number = 1e-12, verbose::Bool = false)
    vars = ["Y"; ["X$i" for i in 1:k]]
    @show(vars)
    edges = [["X$i", "Y"] for i in 1:k]
    push!(edges, ["X$i" for i in 1:k])

    H = Hypergraph(vars, edges)

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
    Y = N - 1
    for y = 0:n-1
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

    Xs = zip(H, ["X$i" for i in 1:k])
    XsY = zip(H, [["X$i" for i in 1:k]; "Y"])
    Y = zip(H, ["Y"])

    Δ = isodd(k) ?
        1 - 1 / ((ω - 2) * floor((k - 1)/2) + ceil((k - 1)/2) + 1) :
        1 - 1 / ((ω - 1) * (k - 1)/2 + 1)
    @show(Δ)

    @constraint(model, -ϵ ≤ h[Y] - 1.0 + Δ / (k-1) ≤ ϵ)
    for i = 1:k
        XiY = zip(H, ["X$i", "Y"])
        Xi = zip(H, ["X$i"])
        @constraint(model, -ϵ ≤ h[XiY] - h[Xi] - Δ ≤ ϵ)
        @constraint(model, -ϵ ≤ h[XiY] - 1.0 ≤ ϵ)
        # @constraint(model, -ϵ ≤ h[XiY] - h[Xi] - h[XsY] + h[Xs] ≤ ϵ)
    end

    @constraint(model, -ϵ ≤ h[XsY] - 1.0 - Δ ≤ ϵ)
    @constraint(model, -ϵ ≤ h[Xs] - 1.0 ≤ ϵ)

    # for Z in combinations(["X$i" for i in 1:k], Int(floor(k/2)))
    #     @constraint(model, -ϵ ≤ h[zip(H, Z)] - sum(h[zip(H, [Xi])] for Xi in Z) ≤ ϵ)
    # end

    A = ["X$i" for i in 1:Int(floor(k/2))]
    B = ["X$i" for i in Int(floor(k/2)) + 1:k]
    # @constraint(model, -ϵ ≤ h[zip(H, A)] - sum(h[zip(H, [Xi])] for Xi in A) ≤ ϵ)
    # @constraint(model, -ϵ ≤ h[zip(H, B)] - sum(h[zip(H, [Xi])] for Xi in B) ≤ ϵ)

    # A_size = Int(floor((k-1)/2))

    # for i in 1:k
    #     X1 = "X$i"
    #     others = ["X$j" for j in 1:k if j != i]
    #     for A in combinations(others, A_size)
    #         B = setdiff(others, A)
    #         AX1 = [A; X1]
    #         BX1 = [B; X1]

    #         @constraint(model, -ϵ ≤ h[zip(H, AX1)] - sum(h[zip(H, [Xi])] for Xi in AX1) ≤ ϵ)
    #         @constraint(model, -ϵ ≤ h[zip(H, BX1)] - sum(h[zip(H, [Xi])] for Xi in BX1) ≤ ϵ)
    #         # @constraint(model, -ϵ ≤ h[zip(H, A)] - h[zip(H, AX1)] + h[zip(H, [X1])] ≤ ϵ)
    #         # @constraint(model, -ϵ ≤ h[zip(H, B)] - h[zip(H, BX1)] + h[zip(H, [X1])] ≤ ϵ)
    #         @constraint(model,
    #             -ϵ ≤ (ω - 2) * h[zip(H, A)] + h[zip(H, B)] + h[zip(H, [X1, "Y"])] - 1.0 - Δ ≤ ϵ)
    #     end
    # end

    optimize!(model)
    @assert termination_status(model) == MathOptInterface.OPTIMAL

    polymatroid = Dict{Set{String}, Float64}()
    verbose && println("\nOptimal Primal Solution:")
    sol = value.(h)
    for i = 0:N-1
        polymatroid[unzip(H, i)] = sol[i]
        # verbose && println("$(f(i)) = $(sol[i])")
    end
    println(polymatroid)
    expr = min_elimination_cost(H, ω, verbose)
    println(expr)
    v = eval(expr, polymatroid)
    for arg in expr.args
        v2 = eval(arg, polymatroid)
        if abs(v - v2) < 1e-6
            println()
            println(arg)
            println()
        end
    end
    @assert abs(1 + Δ - v) < 1e-6 """
     - 1+Δ: $(1 + Δ)
     - v: $v
    """

    return polymatroid
end

end
