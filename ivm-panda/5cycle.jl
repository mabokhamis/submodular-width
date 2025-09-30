using ..HypergraphWidths
io = open("5cycle.log", "a")
println(io, repeat("#", 80))
################################################################################
# Extension 1/6:
H = Hypergraph(
    [:A, :B, :C, :D, :E, :Z2, :Z3, :Z4],
    [
        [:A, :B],
        [:B, :C, :Z2, :Z3, :Z4],
        [:C, :D, :Z2],
        [:D, :E, :Z2, :Z3],
        [:A, :E, :Z2, :Z3, :Z4],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "1/6: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 2/6:
H = Hypergraph(
    [:A, :B, :C, :D, :E, :Z2, :Z3, :Z4],
    [
        [:A, :B],
        [:B, :C, :Z2, :Z3, :Z4],
        [:C, :D, :Z2],
        [:D, :E, :Z2, :Z3, :Z4],
        [:A, :E, :Z2, :Z3],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "2/6: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 3/6:
H = Hypergraph(
    [:A, :B, :C, :D, :E, :Z2, :Z3, :Z4],
    [
        [:A, :B],
        [:B, :C, :Z2],
        [:C, :D, :Z2, :Z3, :Z4],
        [:D, :E, :Z2, :Z3, :Z4],
        [:A, :E, :Z2, :Z3],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "3/6: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 4/6:
H = Hypergraph(
    [:A, :B, :C, :D, :E, :Z2, :Z3, :Z4],
    [
        [:A, :B],
        [:B, :C, :Z2],
        [:C, :D, :Z2, :Z3],
        [:D, :E, :Z2, :Z3, :Z4],
        [:A, :E, :Z2, :Z3, :Z4],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "4/6: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 5/6:
H = Hypergraph(
    [:A, :B, :C, :D, :E, :Z2, :Z3, :Z4],
    [
        [:A, :B],
        [:B, :C, :Z2],
        [:C, :D, :Z2, :Z3, :Z4],
        [:D, :E, :Z2, :Z3],
        [:A, :E, :Z2, :Z3, :Z4],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "5/6: 	")
println(io, submodular_width(H))
flush(io)
################################################################################
# Extension 6/6:
H = Hypergraph(
    [:A, :B, :C, :D, :E, :Z2, :Z3, :Z4],
    [
        [:A, :B],
        [:B, :C, :Z2, :Z3],
        [:C, :D, :Z2],
        [:D, :E, :Z2, :Z3, :Z4],
        [:A, :E, :Z2, :Z3, :Z4],
    ];
    weights = [1.0, 1.0, 1.0, 1.0, 1.0]
)


print(io, "6/6: 	")
println(io, submodular_width(H))
flush(io)
close(io)
