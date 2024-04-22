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
    vars = sort(collect(s); by = x -> (last(x) > 0, length(first(x)), to_string(first(x))))

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

function add_basic_submodularities!(ec::EntropyConstraints, V::Var)
    for Z ∈ powerset(sort(collect(V)))
        W = collect(setdiff(V, Z))
        for (i, x) ∈ enumerate(W), (j, y) ∈ enumerate(W)
            if i < j
                c = submodularity(Set([Z; x]), Set([Z; y]))
                add_constraint!(ec, c)
            end
        end
    end
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
    x = pinv(A) * b
    println(A*x - b)
    for (i, v) ∈ enumerate(x)
        abs(v) < 1e-6 && continue
        println("$v * ($(to_string(ec.constraints[i])))")
    end
end

ec = EntropyConstraints()
add_basic_submodularities!(ec, Set([:A, :B, :C]))
target = Dict(
    Set([:A, :B]) => 1.0,
    Set([:A, :C]) => 1.0,
    Set([:B, :C]) => 1.0,
    Set([:A, :B, :C]) => -2.0,
    Set(Symbol[]) => -1.0,
)

express_constraint(ec, target)

end
