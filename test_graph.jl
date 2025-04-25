include("subw.jl")

using .HypergraphWidths

println(repeat("=", 80))
#= H = Hypergraph(
    [1, 2, 3, 4, 5 ,6, 7, 8, 9, 10, 11, 12],
    [[1, 2], 
    [2, 3],
    [3, 4],
    [4, 1],
    [2, 5],
    [5,6],
    [5,7],
    [6,8],
    [7,8],
    [8,10],
    [10,9],
    [10,12],
    [9,11],
    [12,11],
    [11,4]]
) =#

#H = Hypergraph([1,2,3,4,5,6], [[1,2], [2,3], [3, 4], [4, 5], [5, 6], [6, 1]])
#= 
H = Hypergraph([:A, :B, :C, :D, :E, :F, :G, :H], 
                [[:A, :B], 
                [:A, :C], 
                [:A, :D], 
                [:B, :C], 
                [:B, :D], 
                [:C, :D], 
                [:C, :E], 
                [:C, :F], 
                [:D, :E], 
                [:D, :F], 
                [:E, :F], 
                [:E, :G], 
                [:E, :H], 
                [:F, :G], 
                [:F, :H], 
                [:G, :H], 
                 ]) =#

H = Hypergraph([:A, :B, :C, :D, :E, :F, :G, :H], 
                 [[:A, :B], 
                 [:A, :C], 
                 [:A, :D], 
                 [:B, :C], 
                 [:B, :D], 
                 [:B, :E], 
                 [:C, :D], 
                 [:C, :E], 
                 [:D, :E], 
                 [:D, :G], 
                 [:D, :F], 
                 [:E, :F], 
                 [:E, :G], 
                 [:E, :H], 
                 [:F, :G], 
                 [:F, :H], 
                 [:G, :H], 
                  ])
 

#using .HypergraphWidths: split_context

#ctx = [:F]

#println(split_context(H, PseudoTree{Symbol}(), ctx, 1))

@show(H)
fhtd = fractional_hypertree_depth(H)
@show(fhtd)

@show(H)
fhtwc_1 = fractional_hypertree_depth_with_caching(H, 1)
@show(fhtwc_1)

#= @show(H)
fhtwc_2 = fractional_hypertree_depth_with_caching(H, 2)
@show(fhtwc_2)
 =#
# fhtwc_1 = 2.5