include("subw.jl")

using .HypergraphWidths

println(repeat("=", 80))

println("Theorem 6.5: PTCR ≱ TD[PT]")

println("Figure 4:")
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
                 [:C, :D], 
                 [:D, :E], 
                 [:E, :F], 
])

@show(H)
pt_ex_4 = pt_time(H, 1)
@show(pt_ex_4)

@show(H)
ptcr_1_ex_4 = ptcr_time(H, 1)
@show(ptcr_1_ex_4)

println("Figure 5:")
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
                 [:D, :H], 
                 [:E, :F], 
                 [:E, :G], 
                 [:E, :H], 
                 [:F, :G], 
                 [:F, :H], 
                 [:G, :H], 
                 [:G, :I], 
                 [:G, :J], 
                 [:G, :K], 
                 [:H, :I], 
                 [:H, :K],  
                 [:I, :J], 
                 [:I, :K], 
                 [:J, :K], 
])


@show(H)
pt_ex_5 = pt_time(H, 1)
@show(pt_ex_5)

ptcr_1_ex_5 = ptcr_time(H, 1)
@show(ptcr_1_ex_5)

ptcr_1_5_ex_5 = ptcr_time(H, 1.5)
@show(ptcr_1_5_ex_5)

println(repeat("=", 80))

# Theorem 7.4: RPT < PTCR 
H = Hypergraph([:A, :B, :C, :D, :E, :F, :G, :H, :I, :J, :K, :L], 
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
                 [:D, :H], 
                 [:E, :F], 
                 [:E, :G], 
                 [:E, :H], 
                 [:F, :G], 
                 [:F, :H], 
                 [:G, :H], 
                 [:G, :I], 
                 [:G, :J], 
                 [:G, :K], 
                 [:H, :I], 
                 [:H, :K],  
                 [:I, :J], 
                 [:I, :K], 
                 [:J, :K], 
                 [:L, :A], 
                 [:L, :B], 
                 [:L, :C], 
                 [:L, :D], 
                 [:L, :E], 
                 [:L, :F], 
                 [:L, :G], 
                 [:L, :H], 
                 [:L, :I], 
                 [:L, :J], 
                 [:L, :K], 
])

@show(H)
pt_ex_6 = pt_time(H, 1)
@show(pt_ex_6)

ptcr_1_ex_6 = ptcr_time(H, 1)
@show(ptcr_1_ex_6)

ptcr_1_5_ex_6 = ptcr_time(H, 1.5)
@show(ptcr_1_5_ex_6)

using CSV

CSV.write("figure_examples_runtimes.csv", (Algorithm=["PT", "PTCR", "PT", "PTCR","PTCR","PT","PTCR","PTCR"],
                                            Example=["Figure 4", "Figure 4", "Figure 5", "Figure 5", "Figure 5", "Figure 6", "Figure 6", "Figure 6"],
                                            TimeExponent=[pt_ex_4, ptcr_1_ex_4, pt_ex_5, ptcr_1_ex_5, ptcr_1_5_ex_5, pt_ex_5, ptcr_1_ex_6, ptcr_1_5_ex_6],
                                            SpaceExponent=[0.0, 1.0, 0.0, 1.0, 1.5, 0.0, 1.0, 1.5]))

