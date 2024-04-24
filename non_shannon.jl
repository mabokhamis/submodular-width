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

function to_string(Z::Var)
    return "h($(join(sort(collect(Z)))))"
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
    return Term(
        monotonicity,
        Sum(
            X ∪ Y => 1.0,
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

# -------------------------------------------

# ec = EntropyConstraints()
# V = Set([:A, :B, :C])
# add_basic_shannon!(ec, V)

# s = Sum()
# add_h!(s, 1.0, Set([:A, :B]))
# add_h!(s, 1.0, Set([:B, :C]))
# add_h!(s, 1.0, Set([:A, :C]))
# add_h!(s, -2.0, Set([:A, :B, :C]))

# express_constraint(ec, s)

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
express_constraint(ec, s)

end
