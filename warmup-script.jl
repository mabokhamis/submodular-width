import Pkg
Pkg.add("JuMP")
Pkg.add("Clp")
Pkg.add("MathOptInterface")
Pkg.add("Combinatorics")
Pkg.add("DataStructures")

include("subw.jl")
include("isomorphism.jl")
include("multivariate-extensions.jl")
