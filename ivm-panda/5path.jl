using ..HypergraphWidths
io = open("5path.log", "a")
println(io, repeat("#", 80))
################################################################################
# Extension 1/30:
H = Hypergraph(
    [:A, :B, :C, :D, :E, :F, :Z2, :Z3, :Z4],
    [
        [:A, :B],
        [:B, :C, :Z2],
        [:C, :D, :Z2, :Z3, :Z4],
        [:D, :E, :Z2, :Z3],
        [:E, :F, :Z2, :Z3, :Z4],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "1/30: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 2/30:
H = Hypergraph(
    [:A, :B, :C, :D, :E, :F, :Z2, :Z3, :Z4],
    [
        [:A, :B, :Z2, :Z3],
        [:B, :C, :Z2],
        [:C, :D],
        [:D, :E, :Z2, :Z3, :Z4],
        [:E, :F, :Z2, :Z3, :Z4],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "2/30: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 3/30:
H = Hypergraph(
    [:A, :B, :C, :D, :E, :F, :Z2, :Z3, :Z4],
    [
        [:A, :B, :Z2, :Z3, :Z4],
        [:B, :C],
        [:C, :D, :Z2],
        [:D, :E, :Z2, :Z3, :Z4],
        [:E, :F, :Z2, :Z3],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "3/30: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 4/30:
H = Hypergraph(
    [:A, :B, :C, :D, :E, :F, :Z2, :Z3, :Z4],
    [
        [:A, :B, :Z2, :Z3, :Z4],
        [:B, :C, :Z2],
        [:C, :D],
        [:D, :E, :Z2, :Z3, :Z4],
        [:E, :F, :Z2, :Z3],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "4/30: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 5/30:
H = Hypergraph(
    [:A, :B, :C, :D, :E, :F, :Z2, :Z3, :Z4],
    [
        [:A, :B, :Z2],
        [:B, :C, :Z2, :Z3, :Z4],
        [:C, :D],
        [:D, :E, :Z2, :Z3],
        [:E, :F, :Z2, :Z3, :Z4],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "5/30: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 6/30:
H = Hypergraph(
    [:A, :B, :C, :D, :E, :F, :Z2, :Z3, :Z4],
    [
        [:A, :B, :Z2, :Z3, :Z4],
        [:B, :C],
        [:C, :D, :Z2, :Z3],
        [:D, :E, :Z2, :Z3, :Z4],
        [:E, :F, :Z2],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "6/30: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 7/30:
H = Hypergraph(
    [:A, :B, :C, :D, :E, :F, :Z2, :Z3, :Z4],
    [
        [:A, :B],
        [:B, :C, :Z2, :Z3, :Z4],
        [:C, :D, :Z2, :Z3],
        [:D, :E, :Z2, :Z3, :Z4],
        [:E, :F, :Z2],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "7/30: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 8/30:
H = Hypergraph(
    [:A, :B, :C, :D, :E, :F, :Z2, :Z3, :Z4],
    [
        [:A, :B, :Z2, :Z3],
        [:B, :C],
        [:C, :D, :Z2, :Z3, :Z4],
        [:D, :E, :Z2],
        [:E, :F, :Z2, :Z3, :Z4],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "8/30: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 9/30:
H = Hypergraph(
    [:A, :B, :C, :D, :E, :F, :Z2, :Z3, :Z4],
    [
        [:A, :B, :Z2, :Z3, :Z4],
        [:B, :C],
        [:C, :D, :Z2, :Z3, :Z4],
        [:D, :E, :Z2, :Z3],
        [:E, :F, :Z2],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "9/30: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 10/30:
H = Hypergraph(
    [:A, :B, :C, :D, :E, :F, :Z2, :Z3, :Z4],
    [
        [:A, :B, :Z2, :Z3],
        [:B, :C],
        [:C, :D, :Z2],
        [:D, :E, :Z2, :Z3, :Z4],
        [:E, :F, :Z2, :Z3, :Z4],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "10/30: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 11/30:
H = Hypergraph(
    [:A, :B, :C, :D, :E, :F, :Z2, :Z3, :Z4],
    [
        [:A, :B],
        [:B, :C, :Z2],
        [:C, :D, :Z2, :Z3],
        [:D, :E, :Z2, :Z3, :Z4],
        [:E, :F, :Z2, :Z3, :Z4],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "11/30: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 12/30:
H = Hypergraph(
    [:A, :B, :C, :D, :E, :F, :Z2, :Z3, :Z4],
    [
        [:A, :B],
        [:B, :C, :Z2, :Z3, :Z4],
        [:C, :D, :Z2, :Z3, :Z4],
        [:D, :E, :Z2],
        [:E, :F, :Z2, :Z3],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "12/30: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 13/30:
H = Hypergraph(
    [:A, :B, :C, :D, :E, :F, :Z2, :Z3, :Z4],
    [
        [:A, :B],
        [:B, :C, :Z2, :Z3, :Z4],
        [:C, :D, :Z2],
        [:D, :E, :Z2, :Z3],
        [:E, :F, :Z2, :Z3, :Z4],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "13/30: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 14/30:
H = Hypergraph(
    [:A, :B, :C, :D, :E, :F, :Z2, :Z3, :Z4],
    [
        [:A, :B],
        [:B, :C, :Z2, :Z3, :Z4],
        [:C, :D, :Z2, :Z3],
        [:D, :E, :Z2],
        [:E, :F, :Z2, :Z3, :Z4],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "14/30: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 15/30:
H = Hypergraph(
    [:A, :B, :C, :D, :E, :F, :Z2, :Z3, :Z4],
    [
        [:A, :B],
        [:B, :C, :Z2, :Z3],
        [:C, :D, :Z2, :Z3, :Z4],
        [:D, :E, :Z2],
        [:E, :F, :Z2, :Z3, :Z4],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "15/30: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 16/30:
H = Hypergraph(
    [:A, :B, :C, :D, :E, :F, :Z2, :Z3, :Z4],
    [
        [:A, :B],
        [:B, :C, :Z2, :Z3, :Z4],
        [:C, :D, :Z2],
        [:D, :E, :Z2, :Z3, :Z4],
        [:E, :F, :Z2, :Z3],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "16/30: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 17/30:
H = Hypergraph(
    [:A, :B, :C, :D, :E, :F, :Z2, :Z3, :Z4],
    [
        [:A, :B, :Z2],
        [:B, :C],
        [:C, :D, :Z2, :Z3],
        [:D, :E, :Z2, :Z3, :Z4],
        [:E, :F, :Z2, :Z3, :Z4],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "17/30: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 18/30:
H = Hypergraph(
    [:A, :B, :C, :D, :E, :F, :Z2, :Z3, :Z4],
    [
        [:A, :B],
        [:B, :C, :Z2, :Z3],
        [:C, :D, :Z2, :Z3, :Z4],
        [:D, :E, :Z2, :Z3, :Z4],
        [:E, :F, :Z2],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "18/30: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 19/30:
H = Hypergraph(
    [:A, :B, :C, :D, :E, :F, :Z2, :Z3, :Z4],
    [
        [:A, :B],
        [:B, :C, :Z2, :Z3, :Z4],
        [:C, :D, :Z2, :Z3, :Z4],
        [:D, :E, :Z2, :Z3],
        [:E, :F, :Z2],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "19/30: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 20/30:
H = Hypergraph(
    [:A, :B, :C, :D, :E, :F, :Z2, :Z3, :Z4],
    [
        [:A, :B, :Z2],
        [:B, :C, :Z2, :Z3, :Z4],
        [:C, :D],
        [:D, :E, :Z2, :Z3, :Z4],
        [:E, :F, :Z2, :Z3],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "20/30: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 21/30:
H = Hypergraph(
    [:A, :B, :C, :D, :E, :F, :Z2, :Z3, :Z4],
    [
        [:A, :B, :Z2],
        [:B, :C, :Z2, :Z3],
        [:C, :D],
        [:D, :E, :Z2, :Z3, :Z4],
        [:E, :F, :Z2, :Z3, :Z4],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "21/30: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 22/30:
H = Hypergraph(
    [:A, :B, :C, :D, :E, :F, :Z2, :Z3, :Z4],
    [
        [:A, :B, :Z2, :Z3, :Z4],
        [:B, :C],
        [:C, :D, :Z2],
        [:D, :E, :Z2, :Z3],
        [:E, :F, :Z2, :Z3, :Z4],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "22/30: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 23/30:
H = Hypergraph(
    [:A, :B, :C, :D, :E, :F, :Z2, :Z3, :Z4],
    [
        [:A, :B],
        [:B, :C, :Z2],
        [:C, :D, :Z2, :Z3, :Z4],
        [:D, :E, :Z2, :Z3, :Z4],
        [:E, :F, :Z2, :Z3],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "23/30: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 24/30:
H = Hypergraph(
    [:A, :B, :C, :D, :E, :F, :Z2, :Z3, :Z4],
    [
        [:A, :B, :Z2, :Z3, :Z4],
        [:B, :C],
        [:C, :D, :Z2, :Z3],
        [:D, :E, :Z2],
        [:E, :F, :Z2, :Z3, :Z4],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "24/30: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 25/30:
H = Hypergraph(
    [:A, :B, :C, :D, :E, :F, :Z2, :Z3, :Z4],
    [
        [:A, :B],
        [:B, :C, :Z2, :Z3],
        [:C, :D, :Z2],
        [:D, :E, :Z2, :Z3, :Z4],
        [:E, :F, :Z2, :Z3, :Z4],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "25/30: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 26/30:
H = Hypergraph(
    [:A, :B, :C, :D, :E, :F, :Z2, :Z3, :Z4],
    [
        [:A, :B, :Z2],
        [:B, :C],
        [:C, :D, :Z2, :Z3, :Z4],
        [:D, :E, :Z2, :Z3, :Z4],
        [:E, :F, :Z2, :Z3],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "26/30: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 27/30:
H = Hypergraph(
    [:A, :B, :C, :D, :E, :F, :Z2, :Z3, :Z4],
    [
        [:A, :B, :Z2, :Z3, :Z4],
        [:B, :C],
        [:C, :D, :Z2, :Z3, :Z4],
        [:D, :E, :Z2],
        [:E, :F, :Z2, :Z3],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "27/30: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 28/30:
H = Hypergraph(
    [:A, :B, :C, :D, :E, :F, :Z2, :Z3, :Z4],
    [
        [:A, :B, :Z2, :Z3, :Z4],
        [:B, :C, :Z2],
        [:C, :D],
        [:D, :E, :Z2, :Z3],
        [:E, :F, :Z2, :Z3, :Z4],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "28/30: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 29/30:
H = Hypergraph(
    [:A, :B, :C, :D, :E, :F, :Z2, :Z3, :Z4],
    [
        [:A, :B, :Z2],
        [:B, :C],
        [:C, :D, :Z2, :Z3, :Z4],
        [:D, :E, :Z2, :Z3],
        [:E, :F, :Z2, :Z3, :Z4],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "29/30: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 30/30:
H = Hypergraph(
    [:A, :B, :C, :D, :E, :F, :Z2, :Z3, :Z4],
    [
        [:A, :B, :Z2, :Z3],
        [:B, :C],
        [:C, :D, :Z2, :Z3, :Z4],
        [:D, :E, :Z2, :Z3, :Z4],
        [:E, :F, :Z2],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "30/30: 	")
println(io, submodular_width(H))
flush(io)
close(io)
