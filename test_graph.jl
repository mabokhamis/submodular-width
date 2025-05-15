include("subw.jl")

using .HypergraphWidths

println(repeat("=", 80))
#= 
Ring of Rings =#
H = Hypergraph(
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
)

# Example 23 (Improves Upon TDs)

#=
H = Hypergraph([:A, :B, :C, :D, :E, :F], 
[[:A, :B], 
[:A, :C], 
[:A, :D], 
[:A, :E], 
[:A, :F], 
[:B, :C], 
[:B, :D], 
[:B, :E], 
[:B, :F], 
[:F, :E], 
[:E, :D], 
[:D, :C],
 ])

fhtd = 2.5

Best Pseudotree: PseudoTree{Symbol}(:A, Dict{Symbol, Union{Nothing, Symbol}}(:F => :B, :A => nothing, :D => :E, :B => :D, :E => :A, :C => :B), Dict{Symbol, Set{Symbol}}(:F => Set(), :A => Set([:E]), :D => Set([:B]), :B => Set([:F, :C]), :E => Set([:D]), :C => Set()), Dict(:F => 5, :A => 1, :D => 3, :B => 4, :E => 2, :C => 5), Set([:F, :A, :D, :B, :E, :C]), Dict(:F => [:A, :E, :D, :B], :A => [], :D => [:A, :E], :B => [:A, :E, :D], :E => [:A], :C => [:A, :E, :D, :B]), Dict{Symbol, @NamedTuple{I::Set{Symbol}, S::Set{Symbol}}}(:F => (I = Set([:A]), S = Set([:B, :E])), :A => (I = Set(), S = Set()), :D => (I = Set(), S = Set([:A, :E])), :B => (I = Set([:A]), S = Set([:D, :E])), :E => (I = Set(), S = Set([:A])), :C => (I = Set([:A]), S = Set([:D, :B]))))
Best Var Cover: Dict{Any, Any}(:F => :A, :A => nothing, :D => nothing, :B => :A, :E => nothing, :C => :A)
Fractional Hypertree Depth With Caching: Dict{Any, Any}(:F => :A, :A => nothing, :D => nothing, :B => :A, :E => nothing, :C => :A)
fhtwc_1 = 2.0

BEST TD: Set{Symbol}[Set([:F, :A, :D, :B, :E, :C])]
fhtdw_1 = 2.5


BEST TD: Set{Symbol}[Set([:F, :A, :D, :B, :E]), Set([:A, :D, :B, :C])]
fhtdw_1_5 = 2.0

=#

# Example 24 (Interesting Variable Cover)
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
                 ])
                  
fhtd = 3.0

Best Pseudotree: PseudoTree{Symbol}(:G, Dict{Symbol, Union{Nothing, Symbol}}(:F => :G, :H => :F, :G => nothing, :D => :C, :A => :B, :B => :D, :E => :H, :C => :E), Dict{Symbol, Set{Symbol}}(:F => Set([:H]), :H => Set([:E]), :G => Set([:F]), :D => Set([:B]), :A => Set(), :B => Set([:A]), :E => Set([:C]), :C => Set([:D])), Dict(:F => 2, :H => 3, :G => 1, :D => 6, :A => 8, :B => 7, :E => 4, :C => 5), Set([:F, :H, :G, :D, :A, :B, :E, :C]), Dict(:F => [:G], :H => [:G, :F], :G => [], :D => [:G, :F, :H, :E, :C], :A => [:G, :F, :H, :E, :C, :D, :B], :B => [:G, :F, :H, :E, :C, :D], :E => [:G, :F, :H], :C => [:G, :F, :H, :E]), Dict{Symbol, @NamedTuple{I::Set{Symbol}, S::Set{Symbol}}}(:F => (I = Set(), S = Set([:G])), :H => (I = Set(), S = Set([:F, :G])), :G => (I = Set(), S = Set()), :D => (I = Set([:F]), S = Set([:E, :C])), :A => (I = Set([:C]), S = Set([:D, :B])), :B => (I = Set(), S = Set([:D, :C])), :E => (I = Set([:G]), S = Set([:F, :H])), :C => (I = Set(), S = Set([:F, :E]))))
Best Var Cover: Dict{Any, Any}(:F => nothing, :H => nothing, :G => nothing, :D => :C, :A => :B, :B => nothing, :E => :F, :C => nothing)
Fractional Hypertree Depth With Caching: Dict{Any, Any}(:F => nothing, :H => nothing, :G => nothing, :D => :C, :A => :B, :B => nothing, :E => :F, :C => nothing)
fhtwc_1 = 2.0

BEST TD: Set{Symbol}[Set([:F, :H, :G, :E]), Set([:F, :A, :D, :B, :E, :C])]
fhtdw_1 = 2.0
=#


#= 
# Example 26
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
                 [:D, :F], 
                 [:D, :G], 
                 [:E, :F], 
                 [:E, :G], 
                 [:E, :H], 
                 [:F, :G], 
                 [:F, :H], 
                 [:G, :H], 
                  ])
# Result:
fhtd = 2.5

Best Pseudotree: PseudoTree{Symbol}(:F, Dict{Symbol, Union{Nothing, Symbol}}(:F => nothing, :H => :E, :A => :B, :D => :F, :G => :H, :B => :C, :E => :D, :C => :E), Dict{Symbol, Set{Symbol}}(:F => Set([:D]), :H => Set([:G]), :A => Set(), :D => Set([:E]), :G => Set(), :B => Set([:A]), :E => Set([:H, :C]), :C => Set([:B])), Dict(:F => 1, :H => 4, :A => 6, :D => 2, :G => 5, :B => 5, :E => 3, :C => 4), Set([:F, :H, :A, :D, :G, :B, :E, :C]), Dict(:F => [], :H => [:F, :D, :E], :A => [:F, :D, :E, :C, :B], :D => [:F], :G => [:F, :D, :E, :H], :B => [:F, :D, :E, :C], :E => [:F, :D], :C => [:F, :D, :E]), Dict{Symbol, @NamedTuple{I::Set{Symbol}, S::Set{Symbol}}}(:F => (I = Set(), S = Set()), :H => (I = Set([:F]), S = Set([:D, :E])), :A => (I = Set([:D]), S = Set([:B, :C])), :D => (I = Set(), S = Set([:F])), :G => (I = Set([:F, :D]), S = Set([:H, :E])), :B => (I = Set([:D]), S = Set([:E, :C])), :E => (I = Set(), S = Set([:F, :D])), :C => (I = Set(), S = Set([:D, :E]))))
Best Var Cover: Dict{Any, Any}(:F => nothing, :H => :F, :A => :D, :D => nothing, :G => :H, :B => :D, :E => nothing, :C => nothing)
Fractional Hypertree Depth With Caching: Dict{Any, Any}(:F => nothing, :H => :F, :A => :D, :D => nothing, :G => :H, :B => :D, :E => nothing, :C => nothing)
fhtwc_1 = 2.5


Best Pseudotree: PseudoTree{Symbol}(:F, Dict{Symbol, Union{Nothing, Symbol}}(:F => nothing, :A => :B, :G => :E, :D => :F, :H => :G, :B => :C, :E => :D, :C => :E), Dict{Symbol, Set{Symbol}}(:F => Set([:D]), :A => Set(), :G => Set([:H]), :D => Set([:E]), :H => Set(), :B => Set([:A]), :E => Set([:G, :C]), :C => Set([:B])), Dict(:F => 1, :A => 6, :G => 4, :D => 2, :H => 5, :B => 5, :E => 3, :C => 4), Set([:F, :A, :G, :D, :H, :B, :E, :C]), Dict(:F => [], :A => [:F, :D, :E, :C, :B], :G => [:F, :D, :E], :D => [:F], :H => [:F, :D, :E, :G], :B => [:F, :D, :E, :C], :E => [:F, :D], :C => [:F, :D, :E]), Dict{Symbol, @NamedTuple{I::Set{Symbol}, S::Set{Symbol}}}(:F => (I = Set(), S = Set()), :A => (I = Set(), S = Set([:D, :B, :C])), :G => (I = Set(), S = Set([:F, :D, :E])), :D => (I = Set(), S = Set([:F])), :H => (I = Set(), S = Set([:F, :G, :E])), :B => (I = Set(), S = Set([:D, :E, :C])), :E => (I = Set(), S = Set([:F, :D])), :C => (I = Set(), S = Set([:D, :E]))))
Best Var Cover: Dict{Any, Any}(:F => nothing, :A => nothing, :G => nothing, :D => nothing, :H => nothing, :B => nothing, :E => nothing, :C => nothing)
Fractional Hypertree Depth With Caching: Dict{Any, Any}(:F => nothing, :A => nothing, :G => nothing, :D => nothing, :H => nothing, :B => nothing, :E => nothing, :C => nothing)
fhtwc_1.5 = 2.0

BEST TD: Set{Symbol}[Set([:A, :D, :B, :E, :C]), Set([:F, :H, :G, :D, :E])]
fhtdw_1 = 2.0
=#

# To Camillo: You can try out different query graphs by changing this function.
H = Hypergraph([:A, :B, :C, :D, :E, :F, :G, :H, :I, :J, :K], 
                 [[:A, :B], 
                 [:A, :C], 
                 [:A, :D], 
                 [:B, :C], 
                 [:B, :D], 
                 [:B, :E], 
                 [:C, :D], 
                 [:C, :E], 
                 [:D, :E], 
                 [:D, :F], 
                 [:D, :G], 
                 [:E, :F], 
                 [:E, :G], 
                 [:E, :H], 
                 [:F, :G], 
                 [:F, :H], 
                 [:G, :H], 
                 [:G, :I], 
                 [:G, :J], 
                 [:H, :I], 
                 [:H, :J], 
                 [:H, :K],  
                 [:I, :J], 
                 [:I, :K], 
                 [:J, :K], 
                  ])

#= 
H = Hypergraph([:A, :B, :C], 
                [[:A, :B], 
                [:B, :C], 
                ]) =#

# This is constant-space pseudo-tree depth
#@show(H)
#ptd = fractional_hypertree_depth(H)
#@show(ptd)

# This is traditional fhtw with no space constraint. It's a bit slow, so fair warning.
#@show(H)
#fhtw = fractional_hypertree_width(H)
#@show(fhtw)

# This is pseudo-tree depth with caching & resets
@show(H)
fhtwc_1 = fractional_hypertree_depth_with_caching(H, 1)
@show(fhtwc_1)

# This is pseudo-tree depth without caching & resets + an outer hypertree decomposition w/ bounded separators (i.e. HTD[PT]). It can also be slow.
#@show(H)
#fhtdw_1 = fractional_hypertree_depth_width(H, 1)
#@show(fhtdw_1)

#= @show(H)
fhtwc_2 = fractional_hypertree_depth_with_caching(H, 2)
@show(fhtwc_2)
 =#
# fhtwc_1 = 2.5