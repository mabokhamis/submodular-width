# submodular-width
Tools to compute the _submodular width_. See [this paper](https://www.cs.bme.hu/~dmarx/papers/marx-csp-jacm.pdf) and [this](https://arxiv.org/abs/1612.02503).


# Usage

The easiest way to use this script is to call the test_graph script with:
```
julia --project=. --threads 12 test_graph.jl 
```
You can modify the test_graph script to try out different query graphs.
