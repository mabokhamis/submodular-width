include("subw.jl")


# This module builds on HypergraphWidths and defines a series of width measures which 
# respect space constraints. Specifically, pt_time(H) computes the runtime of a constant-space algorithm
# called pseudo-trees. ptcr_time(H, s)
module SpaceConstrainedWidths

using ..HypergraphWidths

export pt_time, ptcr_time


mutable struct PseudoTree{T}
    root::Union{T, Nothing}
    parent::Dict{T, Union{T, Nothing}}
    children::Dict{T,Set{T}}
    depth::Dict{T,Int}
    vars::Set{T}
    ancestors::Dict{T,Vector{T}}
    context::Dict{T, NamedTuple{(:I, :S), Tuple{Set{T}, Set{T}}}}
    PseudoTree{T}() where T = new{T}(nothing, Dict(), Dict(), Dict(), Set(), Dict(), Dict())
end

function Base.isequal(P1::PseudoTree, P2::PseudoTree)
    return P1.children == P2.children
end

function Base.hash(P1::PseudoTree, h::UInt)
    return hash(P1.children, h)
end

function add_leaf!(P::PseudoTree{T}, v::T, p::Union{T, Nothing}) where T
    P.children[v] = Set{T}()
    P.parent[v] = p
    push!(P.vars, v)
    if isnothing(P.root)
        P.root = v
        P.ancestors[v] = T[]
        P.depth[v] = 1
        return
    end
    P.depth[v] = P.depth[p] + 1
    P.ancestors[v] = union(P.ancestors[p], [p])
    push!(P.children[p], v)
end

function Base.copy(v::Symbol)
    return v
end

function Base.copy(P::PseudoTree{T}) where T
    new_P = PseudoTree{T}()
    new_P.root = copy(P.root)
    new_P.parent = copy(P.parent)
    new_P.children = deepcopy(P.children)
    new_P.depth = copy(P.depth)
    new_P.vars = copy(P.vars)
    new_P.ancestors = deepcopy(P.ancestors)
    new_P.context = deepcopy(P.context)
    return new_P
end

function valid_extensions(P::PseudoTree{T}, H::Hypergraph{T}, only_add_neighbors) where T
    exts = Tuple{T, T}[]
    H_minus_P = induced_hyper_subgraph(H, setdiff(H.vars, P.vars))
    for h_var in setdiff(H.vars, P.vars)
        potential_parents = if only_add_neighbors
            coincident_vars = union([H.edges[idx] for idx in H._var_edges[h_var]]...)
            P.vars ∩ coincident_vars
        else
            P.vars
        end
        connected_vars = connected_component(H_minus_P, h_var)
        # It should never be beneficial to add a (p_var, h_var) edge to the PT which does not
        # occur in H.
        for p_var in potential_parents
            new_branch = union(P.ancestors[p_var], [p_var])
            is_valid = true
            # If an edge includes pieces of the connected component for `h_var`, then it
            # cannot also connect to a different branch. Otherwise, the two branches would
            # need to connect later on.
            for edge in H.edges
                if !isempty(∩(edge, connected_vars))
                    if ∩(edge, P.vars) != ∩(edge, new_branch)
                        is_valid = false
                        break
                    end
                end
            end
            if is_valid
                push!(exts, (h_var, p_var))
            end
        end
    end
    return exts
end

function enumerate_pseudotrees(H::Hypergraph{T}, cost_bound, cost_func, only_add_neighbors, verbose) where T
    pseudotrees = Set{PseudoTree{T}}()
    for v in H.vars
        var_tree = PseudoTree{T}()
        add_leaf!(var_tree, v, nothing)
        push!(pseudotrees, var_tree)
    end
    for i in 1:(length(H.vars)-1)
        nthreads = Threads.nthreads()
        thread_next_pseudotrees = [Set{PseudoTree{T}}() for i in 1:nthreads]  # Per-thread minima initialized to PseudoTree{T}()
        Threads.@threads for pt in collect(pseudotrees)
            thread_id = Threads.threadid()  # Get current thread's ID (1-based)
            for v_ext in valid_extensions(pt, H, only_add_neighbors)
                new_pt = copy(pt)
                add_leaf!(new_pt, v_ext[1], v_ext[2])
                if isinf(cost_bound) || cost_func(H, new_pt) < cost_bound
                    push!(thread_next_pseudotrees[thread_id], new_pt)
                end
            end
        end
        pseudotrees = union(thread_next_pseudotrees...)
        println("Detected PseudoTrees = $(length(pseudotrees))")
    end
    return pseudotrees
end

function get_descendants(H::Hypergraph{T}, P::PseudoTree{T}, var::T) where T
    descendants = Set{T}()
    for v in P.vars
        if var ∈ P.ancestors[v]
            push!(descendants, v)
        end
    end
    return descendants
end


"""
    pt_time(H, [verbose])

Compute the optimal pseudo-tree runtime for the query hypergraph 'H'.

"""

function pt_time(H::Hypergraph{T}, P::PseudoTree{T}) where T
    leaves = [v for v in P.vars if length(P.children[v]) == 0]
    branches = [∪(P.ancestors[v], [v]) for v in leaves]
    branch_depths = [fractional_edge_cover(H, collect(branch)) for branch in branches]
    return maximum(branch_depths; init = 0)
end

function pt_time(H::Hypergraph{T}, verbose=true) where T
    parent_pseudotrees = enumerate_pseudotrees(H, Inf, pt_time, true, verbose)
    parent_cost = minimum([pt_time(H, P) for P in parent_pseudotrees]; init = Inf)
    all_pseudotrees = enumerate_pseudotrees(H, parent_cost, pt_time, false, verbose)
    return min(parent_cost, minimum([pt_time(H, P) for P in all_pseudotrees]; init = Inf))
end


function get_context(H::Hypergraph{T}, P::PseudoTree{T}, var::T) where T
    ancestors = P.ancestors[var]
    descendants = union(get_descendants(H, P, var), [var])
    ctx = Set{T}()
    for edge in H.edges
        if !isempty(∩(edge, ancestors)) && !isempty(∩(edge, descendants))
            union!(ctx, ∩(edge, ancestors))
        end
    end
    return sort(collect(ctx), by=(x)->P.depth[x])
end

function split_context(H::Hypergraph{T}, P::PseudoTree{T}, ctx::Vector{T}, s) where T
    valid_S = T[]
    valid_I = ctx
    for i in 1:length(ctx)
        S = ctx[(length(ctx)-i+1):end]
        I = ctx[1:(length(ctx)-i)]
        S_space = fractional_edge_cover(H, S)
        if S_space ≤ s
            valid_S = S
            valid_I = I
        else
            break
        end
    end
    return (I=Set(valid_I), S=Set(valid_S))
end

function set_contexts!(H::Hypergraph{T}, P::PseudoTree{T}, s) where T
    for var in P.vars
        full_ctx = get_context(H, P, var)
        P.context[var] = split_context(H, P, full_ctx, s)
    end
end

function var_cover_to_RV_sets(H::Hypergraph{T}, P::PseudoTree{T}, vc_dict::Dict{T, Union{T, Nothing}}) where T
    RV_dict = Dict{T, Set{T}}()
    for var in P.vars
        RVs = ∪(Set{T}(P.context[var].I), P.context[var].S, [var])
        vc = vc_dict[var]
        while !isnothing(vc)
            RVs = ∪(RVs, P.context[vc].I, P.context[vc].S, [vc])
            vc = vc_dict[vc]
        end
        RV_dict[var] = RVs
    end
    return RV_dict
end


# Assuming contexts have already been set, compute the minimum relevant variables set.
function compute_min_max_RV_width(H::Hypergraph{T}, P::PseudoTree{T}) where T
    var_cover = Dict{Tuple{T,T}, Vector{T}}()
    for depth in 1:max(values(P.depth)...)
        cur_gen = [var for var in P.vars if P.depth[var] == depth]
        for B in cur_gen
            for A in ∪(P.ancestors[B], [B])
                var_cover[(A, B)] = ∩(∪(P.ancestors[A], [A]), ∪(P.context[B].I, P.context[B].S, [B]))
                shared_anc = ∩(∪(P.ancestors[A], [A]), P.context[B].I)
                if !isempty(shared_anc)
                    recent_shared_anc = argmax((v)->P.depth[v], shared_anc)
                    var_cover[(A, B)] = var_cover[(A, B)] ∪ var_cover[(recent_shared_anc, P.parent[B])]
                end
            end
        end
    end

    relevant_vars = Dict{T, Vector{T}}(A=>var_cover[(A, A)] for A in P.vars)
    rv_width = maximum([fractional_edge_cover(H, rvs) for rvs in values(relevant_vars)])
    return rv_width, var_cover
end

"""
    ptcr_time(H, s, [verbose])

Compute the fractional hypertree depth of hypergraph 'H'.

"""
function ptcr_time(H::Hypergraph{T}, P::PseudoTree{T}, s) where T
    set_contexts!(H, P, s)
    min_RV_width, relevant_vars = compute_min_max_RV_width(H, P)
    return min_RV_width, relevant_vars, P
end

function ptcr_time(H::Hypergraph{T}, s, verbose=true) where T

    parent_pseudotrees = enumerate_pseudotrees(H, Inf, (H, P)->ptcr_time(H, P, s)[1], true, verbose)
    println("ENUMERATED PSEUDOTREES: ", length(parent_pseudotrees))
    nthreads = Threads.nthreads()
    println("N Threads: $nthreads")
    println("Starting 1st PseudoTree Search")
    count = 0
    thread_minimum_depths = fill(Inf, nthreads)  # Per-thread minima initialized to Inf
    thread_min_var_covers = fill(Dict(), nthreads)  # Per-thread minima initialized to Dict()
    thread_min_pseudotrees = fill(PseudoTree{T}(), nthreads)  # Per-thread minima initialized to PseudoTree{T}()
    Threads.@threads for P in collect(parent_pseudotrees)
        thread_id = Threads.threadid()  # Get current thread's ID (1-based)
        count += 1
        if count % 1000 == 0
            println(float(count)/length(parent_pseudotrees), "% done")
        end
        # Update this thread's local minimum safely
        min_RV_width, best_var_cover, P = ptcr_time(H, P, s)
        if min_RV_width < thread_minimum_depths[thread_id]
            thread_minimum_depths[thread_id] = min_RV_width
            thread_min_var_covers[thread_id] = best_var_cover
            thread_min_pseudotrees[thread_id] = P
        end
    end

    minimum_depth = Inf
    min_var_cover = nothing
    min_pseudotree = nothing
    for i in eachindex(thread_minimum_depths)
        if thread_minimum_depths[i] < minimum_depth
            minimum_depth = thread_minimum_depths[i]
            min_var_cover = thread_min_var_covers[i]
            min_pseudotree = thread_min_pseudotrees[i]
        end
    end
    println("Parent Cost: $minimum_depth")
    println("Parent Pseudotree: $min_pseudotree")
    println("Parent Var Cover: $min_var_cover")
    parent_cost = minimum_depth


    println("Starting 2nd PseudoTree Search")
    all_pseudotrees = enumerate_pseudotrees(H, parent_cost, (H, P)->ptcr_time(H, P, s)[1], false, verbose)
    println("ENUMERATED PSEUDOTREES: ", length(all_pseudotrees))

    count = 0
    thread_minimum_depths = fill(Inf, nthreads)  # Per-thread minima initialized to Inf
    thread_min_var_covers = fill(Dict(), nthreads)  # Per-thread minima initialized to Dict()
    thread_min_pseudotrees = fill(PseudoTree{T}(), nthreads)  # Per-thread minima initialized to PseudoTree{T}()
    Threads.@threads for P in collect(all_pseudotrees)
        thread_id = Threads.threadid()  # Get current thread's ID (1-based)
        count += 1
        if count % 100 == 0
            println(float(count)/length(all_pseudotrees), "% done")
        end
        # Update this thread's local minimum safely
        min_RV_width, best_var_cover, P = ptcr_time(H, P, s)
        if min_RV_width < thread_minimum_depths[thread_id]
            thread_minimum_depths[thread_id] = min_RV_width
            thread_min_var_covers[thread_id] = best_var_cover
            thread_min_pseudotrees[thread_id] = P
        end
    end

    minimum_depth = Inf
    min_var_cover = nothing
    min_pseudotree = nothing
    for i in eachindex(thread_minimum_depths)
        if thread_minimum_depths[i] < minimum_depth
            minimum_depth = thread_minimum_depths[i]
            min_var_cover = thread_min_var_covers[i]
            min_pseudotree = thread_min_pseudotrees[i]
        end
    end
    println("Best Var Cover: $min_var_cover")
    println("Best PseudoTree: $min_pseudotree")
    return min(parent_cost, minimum_depth)    
end


end