This repository is a fork of [https://github.com/mabokhamis/submodular-width](https://github.com/mabokhamis/submodular-width). The previous repository contains code for computing classical width measures like fractional hypertree width and submodular width. This repository adds additional measures that adapt to asymptotic space constraints.

# Usage

To reproduce the results from the paper, you need to install the julia programming language. Then, you can simply run :
```
julia --project=. --threads 12 figure_examples.jl 
```
To try out different query graphs, you can modify the hypergraphs defined in this script.
