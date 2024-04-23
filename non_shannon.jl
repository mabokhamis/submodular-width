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

mutable struct EntropyConstraints
    vars::Vector{Var}
    var_index::Dict{Var,Int}
    constraints::Vector{Sum}

    function EntropyConstraints()
        return new(Var[], Dict{Var,Int}(), Sum[])
    end
end

function add_var!(ec::EntropyConstraints, Z::Var)
    haskey(ec.var_index, Z) && return false
    push!(ec.vars, Z)
    ec.var_index[Z] = length(ec.vars)
    return true
end

function add_constraint!(ec::EntropyConstraints, sum::Sum)
    for (var, _) ∈ sum
        add_var!(ec, var)
    end
    push!(ec.constraints, sum)
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

function submodularity(X::Var, Y::Var)
    return Sum(
        X => 1.0,
        Y => 1.0,
        X ∩ Y => -1.0,
        X ∪ Y => -1.0,
    )
end

function monotonicity(X::Var, Y::Var)
    return Sum(
        X ∪ Y => 1.0,
        X => -1.0,
    )
end

function strictness()
    return Sum(Set(Symbol[]) => -1.0,)
end

function add_basic_submodularities!(ec::EntropyConstraints, V::Var)
    V = sort(collect(V))
    for Z ∈ powerset(V)
        W = setdiff(V, Z)
        for (i, x) ∈ enumerate(W), (j, y) ∈ enumerate(W)
            if i < j
                add_constraint!(ec, submodularity(Set([Z; x]), Set([Z; y])))
            end
        end
    end
end

function add_basic_monotonicities!(ec::EntropyConstraints, V::Var)
    for x ∈ sort(collect(V))
        add_constraint!(ec, monotonicity(setdiff(V, Set((x,))), Set((x,))))
    end
    add_constraint!(ec, Sum(Set(Symbol[]) => 1.0,))
end

function add_strictness!(ec::EntropyConstraints)
    add_constraint!(ec, strictness())
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
        for (v, a) ∈ c
            i = ec.var_index[v]
            @assert A[i, j] == 0.0
            A[i, j] = a
        end
    end
    b = zeros(n)
    for (v, a) ∈ sum
        i = ec.var_index[v]
        b[i] = a
    end
    x = solve(A, b)
    for (i, v) ∈ enumerate(x)
        abs(v) < 1e-8 && continue
        println("$v * ($(to_string(ec.constraints[i])))")
    end
end

ec = EntropyConstraints()
V = Set([:A, :B, :C])
add_basic_submodularities!(ec, V)
add_basic_monotonicities!(ec, V)
add_strictness!(ec)

for c in ec.constraints
    println(to_string(c))
end

# target = Dict(
#     Set([:A, :B]) => 1.0,
#     Set([:A, :C]) => 1.0,
#     Set([:B, :C]) => 1.0,
#     Set([:A, :B, :C]) => -2.0,
# )

# express_constraint(ec, target)

target = Dict(
    Set(Symbol[:A]) => 1.0,
)

express_constraint(ec, target)

end
