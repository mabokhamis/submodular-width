# import Pkg
# Pkg.add("JuMP")
# Pkg.add("Clp")
# Pkg.add("MathOptInterface")
# Pkg.add("Combinatorics")
# Pkg.add("DataStructures")

module NonShannon

using JuMP, Clp, MathOptInterface
using Combinatorics
using DataStructures
using LinearAlgebra

const Var = Set{Symbol}
const Sum = Dict{Var,Float64}

function add!(s1::Sum, s2::Sum)
    for (v, a) ∈ s2
        s1[v] = get(s1, v, 0.0) + a
        if s1[v] ≈ 0.0
            delete!(s1, v)
        end
    end
end

function subtract!(s1::Sum, s2::Sum)
    for (v, a) ∈ s2
        s1[v] = get(s1, v, 0.0) - a
        if s1[v] ≈ 0.0
            delete!(s1, v)
        end
    end
end

@enum TermType submodularity monotonicity strictness copy_term conditional_independence conditional

struct Term
    type::TermType
    sum::Sum
end

mutable struct EntropyConstraints
    vars::Vector{Var}
    var_index::Dict{Var,Int}
    constraints::Vector{Term}

    function EntropyConstraints()
        return new(Var[], Dict{Var,Int}(), Term[])
    end
end

function add_var!(ec::EntropyConstraints, Z::Var)
    haskey(ec.var_index, Z) && return false
    push!(ec.vars, Z)
    ec.var_index[Z] = length(ec.vars)
    return true
end

function add_constraint!(ec::EntropyConstraints, t::Term)
    for (var, _) ∈ t.sum
        add_var!(ec, var)
    end
    push!(ec.constraints, t)
end

function to_string(Y::Var, X::Var = Var())
    Y = Y ∪ X
    return isempty(X) ?
        "h($(join(sort(collect(Y)))))" :
        "h($(join(sort(collect(Y))))|$(join(sort(collect(X)))))"
end

function to_string(s::Sum)
    vars = sort(collect(s); by = x -> (last(x) < 0, length(first(x)), to_string(first(x))))

    function coef(a)
        return if a == 1.0
            ""
        elseif a == -1.0
            "-"
        elseif isinteger(a)
            string(Int(a))
        else
            string(a)
        end
    end

    return join(["$(coef(a))$(to_string(v))" for (v, a) ∈ vars], " + ")
end

function to_string(t::Term)
    return "[$(t.type)] $(to_string(t.sum))"
end

function create_submodularity(X::Var, Y::Var)
    @assert !(X ⊆ Y) && !(Y ⊆ X)
    return Term(
            submodularity,
            Sum(
                X => 1.0,
                Y => 1.0,
                X ∩ Y => -1.0,
                X ∪ Y => -1.0,
            )
    )
end

function create_monotonicity(X::Var, Y::Var)
    Y = Y ∪ X
    @assert X != Y
    return Term(
        monotonicity,
        Sum(
            Y => 1.0,
            X => -1.0,
        )
    )
end

function create_strictness()
    return Term(
        strictness,
        Sum(Var() => -1.0,)
    )
end

function add_basic_submodularities!(ec::EntropyConstraints, V::Var)
    V = sort(collect(V))
    for Z ∈ powerset(V)
        W = setdiff(V, Z)
        for (i, x) ∈ enumerate(W), (j, y) ∈ enumerate(W)
            if i < j
                add_constraint!(ec, create_submodularity(Set([Z; x]), Set([Z; y])))
            end
        end
    end
end

function add_basic_monotonicities!(
    ec::EntropyConstraints, V::Var; include_strictness::Bool = false
)
    for x ∈ sort(collect(V))
        add_constraint!(ec, create_monotonicity(setdiff(V, Set((x,))), Set((x,))))
    end
    include_strictness && add_constraint!(ec, Term(
        strictness,
        Sum(Var() => 1.0,))
    )
end

function add_strictness!(ec::EntropyConstraints)
    add_constraint!(ec, create_strictness())
end

function add_basic_shannon!(
    ec::EntropyConstraints, V::Var; include_strictness::Bool = false
)
    add_basic_submodularities!(ec, V)
    add_basic_monotonicities!(ec, V; include_strictness)
    include_strictness && add_strictness!(ec)
end

function add_copy_lemma!(
    ec::EntropyConstraints, X::Var, Y1::Vector{Symbol}, Y2::Vector{Symbol}
)
    @assert length(Y1) == length(Y2) == length(Set(Y1)) == length(Set(Y2))
    @assert isempty(X ∩ Y1) && isempty(X ∩ Y2) && isempty(Y1 ∩ Y2)
    f = Dict(y1 => Y2[i] for (i, y1) ∈ enumerate(Y1))
    V1 = sort(collect(X ∪ Y1))
    for W1 ∈ powerset(V1)
        W1 ⊆ X && continue
        W2 = (get(f, x, x) for x ∈ W1)
        add_constraint!(ec, Term(
            copy_term,
            Sum(
                Set(W1) => 1.0,
                Set(W2) => -1.0,
            )
        ))
        add_constraint!(ec, Term(
            copy_term,
            Sum(
                Set(W2) => 1.0,
                Set(W1) => -1.0,
            )
        ))
    end
    add_constraint!(ec, Term(
        conditional_independence,
        Sum(
            Set(Y1 ∪ X) => -1.0,
            Set(Y2 ∪ X) => -1.0,
            Set(X ∪ Y1 ∪ Y2) => 1.0,
            X => 1.0,
        )
    ))
end

function solve(A, b)
    model = Model(Clp.Optimizer)
    set_optimizer_attribute(model, "LogLevel", 0)
    @variable(model, x[1:size(A, 2)] >= 0)
    @constraint(model, A * x >= b)
    @objective(model, Min, sum(A * x))
    optimize!(model)
    obj = objective_value(model)
    if obj - sum(b) > 1e-8
        error("Infeasible")
    end
    return value.(x)
end

function express_constraint(ec::EntropyConstraints, sum::Sum)
    n = length(ec.vars)
    m = length(ec.constraints)
    A = zeros(n, m)
    for (j, c) ∈ enumerate(ec.constraints)
        for (v, a) ∈ c.sum
            i = ec.var_index[v]
            @assert A[i, j] == 0.0
            A[i, j] = a
        end
    end
    b = zeros(n)
    for (v, a) ∈ sum
        i = ec.var_index[v]
        @assert b[i] == 0.0
        b[i] = a
    end
    x = solve(A, b)
    for (i, v) ∈ enumerate(x)
        abs(v) < 1e-8 && continue
        println("$v * ($(to_string(ec.constraints[i])))")
    end
    return x
end

function add_h!(s::Sum, c::Float64, Y::Var, X::Var = Var())
    Y = Y ∪ X
    s[X] = get(s, X, 0.0) - c
    s[Y] = get(s, Y, 0.0) + c
end

function add_I!(s::Sum, c::Float64, X::Var, Y::Var, Z::Var)
    X = Z ∪ X
    Y = Z ∪ Y
    W = X ∪ Y
    s[X] = get(s, X, 0.0) + c
    s[Y] = get(s, Y, 0.0) + c
    s[Z] = get(s, Z, 0.0) - c
    s[W] = get(s, W, 0.0) - c
end

function create_conditional(Y::Var, X::Var = Var())
    Y = Y ∪ X
    return Term(
        conditional,
        X == Y ? Sum() : Sum(X => -1.0, Y => 1.0,)
    )
end

function _get_X_Y(t::Term)
    Xs = collect(v for (v, c) ∈ t.sum if c == -1.0)
    Ys = collect(v for (v, c) ∈ t.sum if c == 1.0)
    if length(Xs) == length(Ys) == 1
        return (only(Xs), only(Ys))
    else
        @assert isempty(Xs) && isempty(Ys)
        return (Var(), Var())
    end
end

function _get_diamond(t::Term, Y::Var)
    sides = collect(v for (v, c) ∈ t.sum if c == 1.0)
    tb = collect(v for (v, c) ∈ t.sum if c == -1.0)
    @assert length(sides) == length(tb) == 2
    @assert Y ∈ sides
    X = only(setdiff(sides, Set((Y,))))
    sort!(tb; by = Z -> length(Z))
    (Z, W) = tb
    @assert Z == X ∩ Y
    @assert W == X ∪ Y
    return (X, Y, Z, W)
end

function get_diamond2(t::Term, W::Var)
    sides = collect(v for (v, c) ∈ t.sum if c == -1.0)
    tb = collect(v for (v, c) ∈ t.sum if c == 1.0)
    @assert length(sides) == length(tb) == 2
    @assert W ∈ tb
    sort!(tb; by = Z -> length(Z))
    W == last(tb) || return nothing
    Z = first(tb)
    (X, Y) = sides
    @assert Z == X ∩ Y
    @assert W == X ∪ Y
    return (X, Y, Z, W)
end

function maybe_apply_proof_step(δ::Term, t::Term)
    @assert δ.type == conditional
    (X, Y) = _get_X_Y(δ)
    isempty(X) || return nothing
    get(t.sum, Y, 0.0) == (t.type == conditional ? -1.0 : 1.0) || return nothing
    if t.type == conditional
        (Y2, Z) = _get_X_Y(t)
        @assert Y2 == Y
        println("$(to_string(Y)) + $(to_string(Z, Y)) --> $(to_string(Z))")
        return Term[
            create_conditional(Z, Var()),
        ]
    elseif t.type == monotonicity || t.type == copy_term
        (Z, Y2) = _get_X_Y(t)
        @assert Y2 == Y
        println("$(to_string(Y)) --> $(to_string(Z))")
        return Term[
            create_conditional(Z, Var()),
        ]
    elseif t.type == submodularity
        (X, Y, Z, W) = _get_diamond(t, Y)
        println("$(to_string(Y)) --> $(to_string(Z)) + $(to_string(Y, Z))")
        println("$(to_string(Y, Z)) --> $(to_string(W, X))")
        return Term[
            create_conditional(Z, Var()),
            create_conditional(W, X),
        ]
    elseif t.type == conditional_independence
        r = get_diamond2(t, Y)
        if isnothing(r)
            @warn "$(to_string(δ)) $(to_string(t))"
            return nothing
        end
        (X, Y, Z, W) = r
        println("$(to_string(W)) --> $(to_string(X)) + $(to_string(Y, Z))")
        return Term[
            create_conditional(X, Var()),
            create_conditional(Y, Z),
        ]
    else
        error("Unsupported term type: $(to_string(t))")
    end
end

function maybe_apply_proof_step(conditionals::Vector{Term}, terms::Vector{Term})
    for (i1, δ1) ∈ enumerate(conditionals), (i2, δ2) ∈ enumerate(conditionals)
        i1 == i2 && continue
        r = maybe_apply_proof_step(δ1, δ2)
        isnothing(r) && continue
        conditionals = [δ for (i, δ) ∈ enumerate(conditionals) if i ∉ (i1, i2)]
        append!(conditionals, r)
        return (true, conditionals, terms)
    end
    terms = sort(terms; by = t -> t.type)
    for (i, δ) ∈ enumerate(conditionals), (j, t) ∈ enumerate(terms)
        r = maybe_apply_proof_step(δ, t)
        isnothing(r) && continue
        conditionals = [δ for (k, δ) ∈ enumerate(conditionals) if k != i]
        append!(conditionals, r)
        terms = [t for (k, t) ∈ enumerate(terms) if k != j]
        return (true, conditionals, terms)
    end
    return (false, conditionals, terms)
end

function generate_proof_sequence(conditionals::Vector{Term}, terms::Vector{Term})
    println(repeat("-", 40))
    while true
        (changed, conditionals, terms) = maybe_apply_proof_step(conditionals, terms)
        changed || break
    end
    if !isempty(terms)
        @warn("Proof failed")
    end
    println(repeat("-", 40))
    for δ ∈ conditionals
        println("$(to_string(δ))")
    end
    println(repeat("-", 40))
    for t ∈ terms
        println("$(to_string(t))")
    end
    println(repeat("-", 40))
    s = Sum()
    for t ∈ conditionals
        add!(s, t.sum)
    end
    for t ∈ terms
        subtract!(s, t.sum)
    end
    println("Final sum: $(to_string(s))")
end

# -------------------------------------------

# ec = EntropyConstraints()
# V = Set([:A, :B, :C])
# add_basic_shannon!(ec, V)

# s = Sum()
# add_h!(s, 1.0, Set([:A, :B]))
# add_h!(s, 1.0, Set([:B, :C]))
# add_h!(s, 1.0, Set([:A, :C]))
# add_h!(s, -2.0, Set([:A, :B, :C]))

# x = express_constraint(ec, s)

# terms = [t for (i, t) ∈ enumerate(ec.constraints) for j ∈ 1:Int(round(x[i]))]

# conditionals = [
#     create_conditional(Set([:A, :B])),
#     create_conditional(Set([:B, :C])),
#     create_conditional(Set([:A, :C])),
# ]

# generate_proof_sequence(conditionals, terms)

# -------------------------------------------

# terms = [Term(conditional_independence, Sum(
#     Set([:A, :B, :C]) => 1.0,
#     Set([:B]) => 1.0,
#     Set([:A, :B]) => -1.0,
#     Set([:B, :C]) => -1.0,
# ))]

# conditionals = [
#     create_conditional(Set([:A, :B, :C])),
# ]

# generate_proof_sequence(conditionals, terms)

# -------------------------------------------

# ec = EntropyConstraints()
# add_basic_submodularities!(ec, Set([:X, :Y, :A, :B, :A2, :B2]))
# add_copy_lemma!(ec, Set(Symbol[:X, :Y]), [:A, :B], [:A2, :B2])
# s = Sum()
# add_I!(s, 2.0, Set([:X]), Set([:Y]), Set([:A]))
# add_I!(s, 1.0, Set([:X]), Set([:Y]), Set([:B]))
# add_I!(s, 1.0, Set([:A]), Set([:B]), Set(Symbol[]))
# add_I!(s, 1.0, Set([:A]), Set([:Y]), Set([:X]))
# add_I!(s, 1.0, Set([:A]), Set([:X]), Set([:Y]))
# add_I!(s, -1.0, Set([:X]), Set([:Y]), Set(Symbol[]))
# express_constraint(ec, s)

# -------------------------------------------

ec = EntropyConstraints()
add_basic_submodularities!(ec, Set([:X, :Y, :C, :A, :B, :A2, :B2]))
add_copy_lemma!(ec, Set(Symbol[:X, :Y]), [:A, :B], [:A2, :B2])
s = Sum()

add_h!(s, -11.0, Set([:A, :B, :X, :Y, :C]))

add_h!(s, 3.0, Set([:X, :Y]))
add_h!(s, 3.0, Set([:A, :X]))
add_h!(s, 3.0, Set([:A, :Y]))
add_h!(s, 1.0, Set([:B, :X]))
add_h!(s, 1.0, Set([:B, :Y]))
add_h!(s, 5.0, Set([:C]))

add_h!(s, 1.0, Set([:A, :B, :X, :Y, :C]), Set([:A, :B]))
add_h!(s, 4.0, Set([:A, :B, :X, :Y, :C]), Set([:A, :X, :Y]))
add_h!(s, 1.0, Set([:A, :B, :X, :Y, :C]), Set([:B, :X, :Y]))

add_h!(s, 1.0, Set([:A, :B, :X, :Y, :C]), Set([:A, :C]))
add_h!(s, 2.0, Set([:A, :B, :X, :Y, :C]), Set([:X, :C]))
add_h!(s, 2.0, Set([:A, :B, :X, :Y, :C]), Set([:Y, :C]))
x = express_constraint(ec, s)

terms = [t for (i, t) ∈ enumerate(ec.constraints) for j ∈ 1:Int(round(x[i]))]

conditionals = [
    [create_conditional(Set([:X, :Y])) for i ∈ 1:3];
    [create_conditional(Set([:A, :X])) for i ∈ 1:3];
    [create_conditional(Set([:A, :Y])) for i ∈ 1:3];
    [create_conditional(Set([:B, :X])) for i ∈ 1:1];
    [create_conditional(Set([:B, :Y])) for i ∈ 1:1];
    [create_conditional(Set([:C])) for i ∈ 1:5];

    [create_conditional(Set([:A, :B, :X, :Y, :C]), Set([:A, :B])) for i ∈ 1:1];
    [create_conditional(Set([:A, :B, :X, :Y, :C]), Set([:A, :X, :Y])) for i ∈ 1:4];
    [create_conditional(Set([:A, :B, :X, :Y, :C]), Set([:B, :X, :Y])) for i ∈ 1:1];

    [create_conditional(Set([:A, :B, :X, :Y, :C]), Set([:A, :C])) for i ∈ 1:1];
    [create_conditional(Set([:A, :B, :X, :Y, :C]), Set([:X, :C])) for i ∈ 1:2];
    [create_conditional(Set([:A, :B, :X, :Y, :C]), Set([:Y, :C])) for i ∈ 1:2];
]

generate_proof_sequence(conditionals, terms)

end
