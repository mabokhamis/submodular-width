using ..HypergraphWidths
io = open("6path.log", "a")
println(io, repeat("#", 80))
################################################################################
# Extension 1/68:
H = Hypergraph(
    [:B, :C, :D, :E, :Z2, :Z3, :Z4, :Z5],
    [
        [:D, :E, :Z2, :Z3, :Z4],
        [:B, :Z2],
        [:C, :D, :Z2, :Z3, :Z4, :Z5],
        [:B, :C],
        [:E, :Z2, :Z3, :Z4, :Z5],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "1/68: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 2/68:
H = Hypergraph(
    [:B, :C, :D, :E, :Z2, :Z3, :Z4, :Z5],
    [
        [:E, :Z2, :Z3, :Z4, :Z5],
        [:B, :C],
        [:B, :Z2, :Z3, :Z4, :Z5],
        [:C, :D, :Z2],
        [:D, :E, :Z2, :Z3],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "2/68: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 3/68:
H = Hypergraph(
    [:C, :D, :E, :F, :Z3, :Z4, :Z5],
    [
        [:F, :Z3, :Z4, :Z5],
        [:C, :Z3, :Z4, :Z5],
        [:D, :E, :Z3],
        [:C, :D],
        [:E, :F, :Z3, :Z4],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "3/68: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 4/68:
H = Hypergraph(
    [:B, :C, :D, :E, :Z2, :Z3, :Z4, :Z5],
    [
        [:D, :E, :Z2, :Z3, :Z4],
        [:E, :Z2, :Z3, :Z4, :Z5],
        [:B, :C],
        [:B, :Z2, :Z3, :Z4, :Z5],
        [:C, :D, :Z2],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "4/68: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 5/68:
H = Hypergraph(
    [:B, :C, :D, :E, :Z2, :Z3, :Z4, :Z5],
    [
        [:E, :Z2, :Z3, :Z4, :Z5],
        [:C, :D, :Z2, :Z3, :Z4, :Z5],
        [:B, :C],
        [:B, :Z2, :Z3],
        [:D, :E, :Z2],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "5/68: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 6/68:
H = Hypergraph(
    [:C, :D, :Z3, :Z4, :Z5],
    [
        [:D, :Z3, :Z4, :Z5],
        [:C, :Z3, :Z4, :Z5],
        [:C, :D],
    ];
    weights = [1.0, 1.0, 1.0]
)


print(io, "6/68: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 7/68:
H = Hypergraph(
    [:B, :C, :D, :E, :F, :Z2, :Z3, :Z4, :Z5],
    [
        [:B, :Z2],
        [:C, :D, :Z2, :Z3, :Z4, :Z5],
        [:B, :C],
        [:D, :E, :Z2, :Z3, :Z4, :Z5],
        [:F, :Z2, :Z3, :Z4],
        [:E, :F, :Z2, :Z3],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "7/68: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 8/68:
H = Hypergraph(
    [:B, :C, :D, :E, :Z2, :Z3, :Z4, :Z5],
    [
        [:E, :Z2, :Z3, :Z4],
        [:C, :D, :Z2, :Z3, :Z4, :Z5],
        [:B, :C],
        [:B, :Z2, :Z3, :Z4, :Z5],
        [:D, :E, :Z2],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "8/68: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 9/68:
H = Hypergraph(
    [:C, :D, :E, :F, :Z3, :Z4, :Z5],
    [
        [:F, :Z3, :Z4, :Z5],
        [:C, :Z3, :Z4, :Z5],
        [:C, :D],
        [:D, :E, :Z3, :Z4],
        [:E, :F, :Z3],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "9/68: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 10/68:
H = Hypergraph(
    [:B, :C, :D, :E, :F, :Z2, :Z3, :Z4, :Z5],
    [
        [:B, :C, :Z2],
        [:D, :E, :Z2, :Z3, :Z4, :Z5],
        [:B, :Z2, :Z3],
        [:E, :F, :Z2, :Z3, :Z4],
        [:F, :Z2, :Z3, :Z4, :Z5],
        [:C, :D],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "10/68: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 11/68:
H = Hypergraph(
    [:C, :D, :E, :Z3, :Z4, :Z5],
    [
        [:E, :Z3, :Z4, :Z5],
        [:C, :Z3, :Z4, :Z5],
        [:C, :D],
        [:D, :E, :Z3, :Z4],
    ];
    weights = [1.0, 1.0, 1.0, 1.0]
)


print(io, "11/68: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 12/68:
H = Hypergraph(
    [:B, :C, :D, :E, :F, :Z2, :Z3, :Z4, :Z5],
    [
        [:F, :Z2, :Z3],
        [:D, :E, :Z2, :Z3, :Z4, :Z5],
        [:B, :Z2, :Z3, :Z4, :Z5],
        [:B, :C, :Z2, :Z3, :Z4],
        [:E, :F, :Z2],
        [:C, :D],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "12/68: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 13/68:
H = Hypergraph(
    [:C, :D, :E, :Z3, :Z4, :Z5],
    [
        [:E, :Z3, :Z4, :Z5],
        [:C, :Z3, :Z4, :Z5],
        [:D, :E, :Z3],
        [:C, :D],
    ];
    weights = [1.0, 1.0, 1.0, 1.0]
)


print(io, "13/68: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 14/68:
H = Hypergraph(
    [:B, :C, :D, :E, :F, :Z2, :Z3, :Z4, :Z5],
    [
        [:B, :Z2, :Z3, :Z4],
        [:D, :E, :Z2, :Z3, :Z4, :Z5],
        [:B, :C],
        [:C, :D, :Z2],
        [:E, :F, :Z2, :Z3],
        [:F, :Z2, :Z3, :Z4, :Z5],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "14/68: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 15/68:
H = Hypergraph(
    [:B, :C, :D, :E, :F, :Z2, :Z3, :Z4, :Z5],
    [
        [:D, :E, :Z2, :Z3, :Z4],
        [:C, :D, :Z2, :Z3, :Z4, :Z5],
        [:B, :C],
        [:B, :Z2, :Z3],
        [:E, :F, :Z2],
        [:F, :Z2, :Z3, :Z4, :Z5],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "15/68: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 16/68:
H = Hypergraph(
    [:B, :C, :D, :E, :Z2, :Z3, :Z4, :Z5],
    [
        [:D, :E, :Z2, :Z3, :Z4],
        [:E, :Z2, :Z3, :Z4, :Z5],
        [:C, :D, :Z2, :Z3, :Z4, :Z5],
        [:B, :C],
        [:B, :Z2, :Z3],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "16/68: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 17/68:
H = Hypergraph(
    [:B, :C, :D, :E, :Z2, :Z3, :Z4, :Z5],
    [
        [:B, :Z2, :Z3, :Z4],
        [:E, :Z2, :Z3, :Z4, :Z5],
        [:C, :D, :Z2, :Z3, :Z4, :Z5],
        [:B, :C],
        [:D, :E, :Z2],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "17/68: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 18/68:
H = Hypergraph(
    [:B, :C, :D, :E, :F, :Z2, :Z3, :Z4, :Z5],
    [
        [:F, :Z2, :Z3],
        [:D, :E, :Z2, :Z3, :Z4],
        [:C, :D, :Z2, :Z3, :Z4, :Z5],
        [:B, :C],
        [:B, :Z2, :Z3, :Z4, :Z5],
        [:E, :F, :Z2],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "18/68: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 19/68:
H = Hypergraph(
    [:B, :C, :D, :E, :Z2, :Z3, :Z4, :Z5],
    [
        [:E, :Z2, :Z3, :Z4, :Z5],
        [:B, :C],
        [:B, :Z2, :Z3, :Z4, :Z5],
        [:C, :D, :Z2, :Z3, :Z4],
        [:D, :E, :Z2, :Z3],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "19/68: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 20/68:
H = Hypergraph(
    [:B, :C, :D, :E, :F, :Z2, :Z3, :Z4, :Z5],
    [
        [:B, :C],
        [:B, :Z2, :Z3, :Z4, :Z5],
        [:C, :D, :Z2],
        [:E, :F, :Z2, :Z3, :Z4],
        [:F, :Z2, :Z3, :Z4, :Z5],
        [:D, :E, :Z2, :Z3],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "20/68: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 21/68:
H = Hypergraph(
    [:B, :C, :D, :E, :F, :Z2, :Z3, :Z4, :Z5],
    [
        [:C, :D, :Z2, :Z3, :Z4, :Z5],
        [:B, :C],
        [:B, :Z2, :Z3, :Z4, :Z5],
        [:F, :Z2, :Z3, :Z4],
        [:D, :E, :Z2],
        [:E, :F, :Z2, :Z3],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "21/68: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 22/68:
H = Hypergraph(
    [:D, :E, :Z4, :Z5],
    [
        [:D, :Z4, :Z5],
        [:D, :E],
        [:E, :Z4, :Z5],
    ];
    weights = [1.0, 1.0, 1.0]
)


print(io, "22/68: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 23/68:
H = Hypergraph(
    [:B, :C, :D, :E, :F, :Z2, :Z3, :Z4, :Z5],
    [
        [:B, :C, :Z2],
        [:B, :Z2, :Z3, :Z4, :Z5],
        [:E, :F, :Z2, :Z3, :Z4],
        [:F, :Z2, :Z3, :Z4, :Z5],
        [:C, :D],
        [:D, :E, :Z2, :Z3],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "23/68: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 24/68:
H = Hypergraph(
    [:B, :C, :Z2, :Z3, :Z4, :Z5],
    [
        [:C, :Z2, :Z3, :Z4, :Z5],
        [:B, :C],
        [:B, :Z2, :Z3, :Z4, :Z5],
    ];
    weights = [1.0, 1.0, 1.0]
)


print(io, "24/68: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 25/68:
H = Hypergraph(
    [:B, :C, :D, :E, :F, :Z2, :Z3, :Z4, :Z5],
    [
        [:D, :E, :Z2, :Z3, :Z4],
        [:B, :Z2],
        [:C, :D, :Z2, :Z3, :Z4, :Z5],
        [:B, :C],
        [:E, :F, :Z2, :Z3],
        [:F, :Z2, :Z3, :Z4, :Z5],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "25/68: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 26/68:
H = Hypergraph(
    [:B, :C, :D, :E, :F, :Z2, :Z3, :Z4, :Z5],
    [
        [:D, :E, :Z2, :Z3, :Z4, :Z5],
        [:B, :C],
        [:C, :D, :Z2],
        [:B, :Z2, :Z3],
        [:E, :F, :Z2, :Z3, :Z4],
        [:F, :Z2, :Z3, :Z4, :Z5],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "26/68: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 27/68:
H = Hypergraph(
    [:B, :C, :D, :E, :Z2, :Z3, :Z4, :Z5],
    [
        [:B, :Z2],
        [:C, :D, :Z2, :Z3, :Z4, :Z5],
        [:B, :C],
        [:E, :Z2, :Z3, :Z4, :Z5],
        [:D, :E, :Z2, :Z3],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "27/68: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 28/68:
H = Hypergraph(
    [:B, :C, :D, :E, :F, :Z2, :Z3, :Z4, :Z5],
    [
        [:B, :Z2, :Z3, :Z4],
        [:C, :D, :Z2, :Z3, :Z4, :Z5],
        [:B, :C],
        [:D, :E, :Z2],
        [:E, :F, :Z2, :Z3],
        [:F, :Z2, :Z3, :Z4, :Z5],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "28/68: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 29/68:
H = Hypergraph(
    [:B, :C, :D, :E, :F, :Z2, :Z3, :Z4, :Z5],
    [
        [:D, :E, :Z2, :Z3, :Z4, :Z5],
        [:B, :Z2, :Z3, :Z4, :Z5],
        [:F, :Z2, :Z3, :Z4],
        [:B, :C, :Z2, :Z3],
        [:E, :F, :Z2],
        [:C, :D],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "29/68: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 30/68:
H = Hypergraph(
    [:B, :C, :D, :E, :Z2, :Z3, :Z4, :Z5],
    [
        [:E, :Z2, :Z3, :Z4, :Z5],
        [:B, :C],
        [:B, :Z2, :Z3, :Z4, :Z5],
        [:D, :E, :Z2],
        [:C, :D, :Z2, :Z3, :Z4],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "30/68: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 31/68:
H = Hypergraph(
    [:D, :E, :F, :Z4, :Z5],
    [
        [:D, :Z4, :Z5],
        [:D, :E],
        [:F, :Z4, :Z5],
        [:E, :F, :Z4],
    ];
    weights = [1.0, 1.0, 1.0, 1.0]
)


print(io, "31/68: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 32/68:
H = Hypergraph(
    [:C, :D, :E, :F, :Z3, :Z4, :Z5],
    [
        [:D, :E, :Z3, :Z4, :Z5],
        [:C, :Z3, :Z4, :Z5],
        [:C, :D],
        [:F, :Z3, :Z4],
        [:E, :F, :Z3],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "32/68: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 33/68:
H = Hypergraph(
    [:B, :C, :D, :E, :F, :Z2, :Z3, :Z4, :Z5],
    [
        [:B, :Z2],
        [:D, :E, :Z2, :Z3, :Z4, :Z5],
        [:B, :C],
        [:E, :F, :Z2, :Z3],
        [:F, :Z2, :Z3, :Z4, :Z5],
        [:C, :D, :Z2, :Z3, :Z4],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "33/68: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 34/68:
H = Hypergraph(
    [:B, :C, :D, :E, :F, :Z2, :Z3, :Z4, :Z5],
    [
        [:B, :C],
        [:B, :Z2, :Z3, :Z4, :Z5],
        [:E, :F, :Z2],
        [:F, :Z2, :Z3, :Z4, :Z5],
        [:C, :D, :Z2, :Z3, :Z4],
        [:D, :E, :Z2, :Z3],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "34/68: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 35/68:
H = Hypergraph(
    [:E, :F, :Z5],
    [
        [:E, :Z5],
        [:F, :Z5],
        [:E, :F],
    ];
    weights = [1.0, 1.0, 1.0]
)


print(io, "35/68: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 36/68:
H = Hypergraph(
    [:B, :C, :D, :E, :F, :Z2, :Z3, :Z4, :Z5],
    [
        [:C, :D, :Z2, :Z3, :Z4, :Z5],
        [:B, :C],
        [:B, :Z2, :Z3],
        [:D, :E, :Z2],
        [:E, :F, :Z2, :Z3, :Z4],
        [:F, :Z2, :Z3, :Z4, :Z5],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "36/68: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 37/68:
H = Hypergraph(
    [:B, :C, :D, :E, :F, :Z2, :Z3, :Z4, :Z5],
    [
        [:B, :C],
        [:B, :Z2, :Z3, :Z4, :Z5],
        [:D, :E, :Z2],
        [:E, :F, :Z2, :Z3],
        [:F, :Z2, :Z3, :Z4, :Z5],
        [:C, :D, :Z2, :Z3, :Z4],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "37/68: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 38/68:
H = Hypergraph(
    [:B, :C, :D, :E, :F, :Z2, :Z3, :Z4, :Z5],
    [
        [:B, :Z2, :Z3, :Z4],
        [:B, :C, :Z2],
        [:D, :E, :Z2, :Z3, :Z4, :Z5],
        [:E, :F, :Z2, :Z3],
        [:F, :Z2, :Z3, :Z4, :Z5],
        [:C, :D],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "38/68: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 39/68:
H = Hypergraph(
    [:C, :D, :E, :F, :Z2, :Z3, :Z4, :Z5],
    [
        [:F, :Z2, :Z3],
        [:C, :Z2, :Z3, :Z4, :Z5],
        [:D, :E, :Z2, :Z3, :Z4, :Z5],
        [:E, :F, :Z2],
        [:C, :D],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "39/68: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 40/68:
H = Hypergraph(
    [:B, :C, :D, :E, :F, :Z2, :Z3, :Z4, :Z5],
    [
        [:B, :Z2, :Z3, :Z4],
        [:F, :Z2, :Z3],
        [:C, :D, :Z2, :Z3, :Z4, :Z5],
        [:B, :C],
        [:D, :E, :Z2, :Z3, :Z4, :Z5],
        [:E, :F, :Z2],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "40/68: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 41/68:
H = Hypergraph(
    [:B, :C, :D, :E, :F, :Z2, :Z3, :Z4, :Z5],
    [
        [:C, :D, :Z2, :Z3, :Z4, :Z5],
        [:B, :C],
        [:B, :Z2, :Z3, :Z4, :Z5],
        [:F, :Z2, :Z3, :Z4],
        [:E, :F, :Z2],
        [:D, :E, :Z2, :Z3],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "41/68: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 42/68:
H = Hypergraph(
    [:C],
    [
        [:C],
    ];
    weights = [1.0]
)


print(io, "42/68: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 43/68:
H = Hypergraph(
    [:B, :C, :D, :E, :Z2, :Z3, :Z4, :Z5],
    [
        [:E, :Z2, :Z3, :Z4],
        [:C, :D, :Z2, :Z3, :Z4, :Z5],
        [:B, :C],
        [:B, :Z2, :Z3, :Z4, :Z5],
        [:D, :E, :Z2, :Z3],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "43/68: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 44/68:
H = Hypergraph(
    [:C, :D, :E, :F, :Z3, :Z4, :Z5],
    [
        [:F, :Z3, :Z4, :Z5],
        [:D, :E, :Z3, :Z4, :Z5],
        [:C, :Z3],
        [:C, :D],
        [:E, :F, :Z3, :Z4],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "44/68: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 45/68:
H = Hypergraph(
    [:B, :C, :D, :E, :F, :Z2, :Z3, :Z4, :Z5],
    [
        [:D, :E, :Z2, :Z3, :Z4],
        [:B, :C],
        [:B, :Z2, :Z3, :Z4, :Z5],
        [:C, :D, :Z2, :Z3],
        [:E, :F, :Z2],
        [:F, :Z2, :Z3, :Z4, :Z5],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "45/68: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 46/68:
H = Hypergraph(
    [:B, :C, :D, :E, :F, :Z2, :Z3, :Z4, :Z5],
    [
        [:B, :C, :Z2],
        [:D, :E, :Z2, :Z3, :Z4, :Z5],
        [:B, :Z2, :Z3, :Z4, :Z5],
        [:F, :Z2, :Z3, :Z4],
        [:E, :F, :Z2, :Z3],
        [:C, :D],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "46/68: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 47/68:
H = Hypergraph(
    [:B, :C, :D, :E, :Z2, :Z3, :Z4, :Z5],
    [
        [:D, :E, :Z2, :Z3, :Z4],
        [:E, :Z2, :Z3, :Z4, :Z5],
        [:B, :C],
        [:B, :Z2, :Z3, :Z4, :Z5],
        [:C, :D, :Z2, :Z3],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "47/68: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 48/68:
H = Hypergraph(
    [:B, :C, :D, :E, :F, :Z2, :Z3, :Z4, :Z5],
    [
        [:D, :E, :Z2, :Z3, :Z4],
        [:B, :C],
        [:B, :Z2, :Z3, :Z4, :Z5],
        [:C, :D, :Z2],
        [:E, :F, :Z2, :Z3],
        [:F, :Z2, :Z3, :Z4, :Z5],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "48/68: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 49/68:
H = Hypergraph(
    [:B, :C, :D, :E, :F, :Z2, :Z3, :Z4, :Z5],
    [
        [:D, :E, :Z2, :Z3, :Z4, :Z5],
        [:B, :C],
        [:B, :Z2, :Z3, :Z4, :Z5],
        [:F, :Z2, :Z3, :Z4],
        [:C, :D, :Z2, :Z3],
        [:E, :F, :Z2],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "49/68: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 50/68:
H = Hypergraph(
    [:B, :C, :D, :E, :F, :Z2, :Z3, :Z4, :Z5],
    [
        [:B, :Z2, :Z3, :Z4],
        [:D, :E, :Z2, :Z3, :Z4, :Z5],
        [:B, :C, :Z2, :Z3],
        [:E, :F, :Z2],
        [:F, :Z2, :Z3, :Z4, :Z5],
        [:C, :D],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "50/68: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 51/68:
H = Hypergraph(
    [:B, :C, :D, :E, :F, :Z2, :Z3, :Z4, :Z5],
    [
        [:C, :D, :Z2, :Z3, :Z4, :Z5],
        [:B, :C],
        [:D, :E, :Z2, :Z3, :Z4, :Z5],
        [:F, :Z2, :Z3, :Z4],
        [:B, :Z2, :Z3],
        [:E, :F, :Z2],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "51/68: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 52/68:
H = Hypergraph(
    [:B, :C, :D, :E, :F, :Z2, :Z3, :Z4, :Z5],
    [
        [:D, :E, :Z2, :Z3, :Z4, :Z5],
        [:B, :C],
        [:B, :Z2, :Z3, :Z4, :Z5],
        [:F, :Z2, :Z3, :Z4],
        [:C, :D, :Z2],
        [:E, :F, :Z2, :Z3],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "52/68: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 53/68:
H = Hypergraph(
    [:B, :C, :D, :E, :F, :Z2, :Z3, :Z4, :Z5],
    [
        [:B, :C],
        [:B, :Z2, :Z3, :Z4, :Z5],
        [:C, :D, :Z2, :Z3],
        [:E, :F, :Z2, :Z3, :Z4],
        [:D, :E, :Z2],
        [:F, :Z2, :Z3, :Z4, :Z5],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "53/68: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 54/68:
H = Hypergraph(
    [:B, :C, :D, :E, :F, :Z2, :Z3, :Z4, :Z5],
    [
        [:B, :C, :Z2],
        [:D, :E, :Z2, :Z3, :Z4],
        [:B, :Z2, :Z3, :Z4, :Z5],
        [:E, :F, :Z2, :Z3],
        [:F, :Z2, :Z3, :Z4, :Z5],
        [:C, :D],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "54/68: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 55/68:
H = Hypergraph(
    [:B, :C, :D, :Z2, :Z3, :Z4, :Z5],
    [
        [:B, :C],
        [:B, :Z2, :Z3, :Z4, :Z5],
        [:C, :D, :Z2, :Z3],
        [:D, :Z2, :Z3, :Z4, :Z5],
    ];
    weights = [1.0, 1.0, 1.0, 1.0]
)


print(io, "55/68: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 56/68:
H = Hypergraph(
    [:B, :C, :D, :Z2, :Z3, :Z4, :Z5],
    [
        [:B, :C],
        [:B, :Z2, :Z3, :Z4, :Z5],
        [:D, :Z2, :Z3, :Z4, :Z5],
        [:C, :D, :Z2, :Z3, :Z4],
    ];
    weights = [1.0, 1.0, 1.0, 1.0]
)


print(io, "56/68: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 57/68:
H = Hypergraph(
    [:C, :D, :E, :F, :Z3, :Z4, :Z5],
    [
        [:F, :Z3, :Z4, :Z5],
        [:D, :E, :Z3, :Z4, :Z5],
        [:C, :D],
        [:C, :Z3, :Z4],
        [:E, :F, :Z3],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "57/68: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 58/68:
H = Hypergraph(
    [:B, :C, :D, :E, :Z2, :Z3, :Z4, :Z5],
    [
        [:E, :Z2, :Z3, :Z4, :Z5],
        [:B, :C],
        [:B, :Z2, :Z3, :Z4, :Z5],
        [:C, :D, :Z2, :Z3],
        [:D, :E, :Z2],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "58/68: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 59/68:
H = Hypergraph(
    [:B, :C, :D, :E, :F, :Z2, :Z3, :Z4, :Z5],
    [
        [:B, :Z2, :Z3, :Z4],
        [:C, :D, :Z2, :Z3, :Z4, :Z5],
        [:B, :C],
        [:E, :F, :Z2],
        [:F, :Z2, :Z3, :Z4, :Z5],
        [:D, :E, :Z2, :Z3],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "59/68: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 60/68:
H = Hypergraph(
    [:B, :C, :D, :E, :F, :Z2, :Z3, :Z4, :Z5],
    [
        [:B, :Z2],
        [:C, :D, :Z2, :Z3, :Z4, :Z5],
        [:B, :C],
        [:E, :F, :Z2, :Z3, :Z4],
        [:F, :Z2, :Z3, :Z4, :Z5],
        [:D, :E, :Z2, :Z3],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "60/68: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 61/68:
H = Hypergraph(
    [:B, :C, :D, :E, :F, :Z2, :Z3, :Z4, :Z5],
    [
        [:D, :E, :Z2, :Z3, :Z4, :Z5],
        [:B, :C],
        [:B, :Z2, :Z3],
        [:E, :F, :Z2],
        [:F, :Z2, :Z3, :Z4, :Z5],
        [:C, :D, :Z2, :Z3, :Z4],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "61/68: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 62/68:
H = Hypergraph(
    [:B, :C, :D, :E, :F, :Z2, :Z3, :Z4, :Z5],
    [
        [:D, :E, :Z2, :Z3, :Z4],
        [:B, :Z2, :Z3, :Z4, :Z5],
        [:B, :C, :Z2, :Z3],
        [:E, :F, :Z2],
        [:F, :Z2, :Z3, :Z4, :Z5],
        [:C, :D],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "62/68: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 63/68:
H = Hypergraph(
    [:B, :C, :D, :E, :Z2, :Z3, :Z4, :Z5],
    [
        [:B, :C, :Z2],
        [:E, :Z2, :Z3, :Z4, :Z5],
        [:B, :Z2, :Z3, :Z4, :Z5],
        [:C, :D],
        [:D, :E, :Z2, :Z3],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "63/68: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 64/68:
H = Hypergraph(
    [:B, :C, :D, :E, :F, :Z2, :Z3, :Z4, :Z5],
    [
        [:B, :Z2, :Z3, :Z4],
        [:D, :E, :Z2, :Z3, :Z4, :Z5],
        [:B, :C],
        [:C, :D, :Z2, :Z3],
        [:E, :F, :Z2],
        [:F, :Z2, :Z3, :Z4, :Z5],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "64/68: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 65/68:
H = Hypergraph(
    [:B, :C, :D, :E, :F, :Z2, :Z3, :Z4, :Z5],
    [
        [:B, :Z2],
        [:D, :E, :Z2, :Z3, :Z4, :Z5],
        [:B, :C],
        [:C, :D, :Z2, :Z3],
        [:E, :F, :Z2, :Z3, :Z4],
        [:F, :Z2, :Z3, :Z4, :Z5],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "65/68: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 66/68:
H = Hypergraph(
    [:B, :C, :D, :Z2, :Z3, :Z4, :Z5],
    [
        [:B, :C],
        [:B, :Z2, :Z3, :Z4, :Z5],
        [:C, :D, :Z2],
        [:D, :Z2, :Z3, :Z4, :Z5],
    ];
    weights = [1.0, 1.0, 1.0, 1.0]
)


print(io, "66/68: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 67/68:
H = Hypergraph(
    [:B, :C, :D, :E, :F, :Z2, :Z3, :Z4, :Z5],
    [
        [:F, :Z2, :Z3],
        [:D, :E, :Z2, :Z3, :Z4, :Z5],
        [:B, :C],
        [:B, :Z2, :Z3, :Z4, :Z5],
        [:E, :F, :Z2],
        [:C, :D, :Z2, :Z3, :Z4],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "67/68: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 68/68:
H = Hypergraph(
    [:B, :C, :D, :E, :Z2, :Z3, :Z4, :Z5],
    [
        [:B, :Z2, :Z3, :Z4],
        [:E, :Z2, :Z3, :Z4, :Z5],
        [:C, :D, :Z2, :Z3, :Z4, :Z5],
        [:B, :C],
        [:D, :E, :Z2, :Z3],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "68/68: 	")
println(io, submodular_width(H))
flush(io)
close(io)
