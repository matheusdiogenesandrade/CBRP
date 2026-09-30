using Logging
using CPLEX

"""Path-CBRP SEC: (S, i, j, k_yi, k_yj) with y_meta[k]=(b,i)."""
const PathSubtourCut = Tuple{Set{Int},Int,Int,Int,Int}

"""Depot-rooted Path-CBRP SEC `sum_{a in δ⁺(S)} x_a >= y_k`: (S, k) with depot ∉ S, y_meta[k][2] ∈ S."""
const PathDepotSubtourCut = Tuple{Set{Int},Int}

"""Canonical user-cut key `(id(S), min(ky_i, ky_j), max(ky_i, ky_j))`; `S` is interned in the stats."""
const PathUserCutKey = Tuple{Int,Int,Int}

"""
Statistics collected by the Path SEC separation callback (user + lazy cuts).

`n_user_cuts` counts cuts actually submitted to CPLEX. `n_user_found` counts every cut returned by
max-flow separation; `n_user_dup` those whose key was already submitted earlier in the solve
(skipped iff `dedup_user_cuts`). `n_user_nonviolated` counts submitted cuts whose violation on the
exact (untruncated) `x` is `<= epsilon`. `user_viol_*` summarise the exact violation of submitted cuts.
`user_cut_max_depth` limits fractional separation (`-1` = unlimited); `n_user_skipped_depth` counts
RELAXATION calls skipped because the node was deeper.
"""
mutable struct PathSubtourCallbackStats
    n_user_cuts::Int
    n_lazy_cuts::Int
    sep_time::Float64
    dedup_user_cuts::Bool
    n_user_found::Int
    n_user_dup::Int
    n_user_nonviolated::Int
    n_user_rounds::Int
    user_viol_min::Float64
    user_viol_max::Float64
    user_viol_sum::Float64
    s_ids::Dict{Vector{Int},Int}
    seen_user_cuts::Set{PathUserCutKey}
    user_cut_max_depth::Int
    n_user_skipped_depth::Int
end

PathSubtourCallbackStats(
    n_user_cuts::Int,
    n_lazy_cuts::Int,
    sep_time::Float64;
    dedup_user_cuts::Bool=false,
    user_cut_max_depth::Int=-1,
) =
    PathSubtourCallbackStats(
        n_user_cuts, n_lazy_cuts, sep_time, dedup_user_cuts,
        0, 0, 0, 0, Inf, -Inf, 0.0,
        Dict{Vector{Int},Int}(), Set{PathUserCutKey}(),
        user_cut_max_depth, 0,
    )

"""
Require callback SEC separation when compact arc MTZ is disabled.
"""
function validatePathCbrpNoMtzSeparation!(app::Dict{String,Any})::Nothing
    get(app, "no-path-cbrp-mtz", false) || return nothing
    sep_mode::String = get(app, "subcycle-separation", "none")
    sep_engine::String = get(app, "subcycle-separation-engine", "root")
    if sep_mode == "none" || sep_engine != "callback"
        throw(ArgumentError(
            "no-path-cbrp-mtz requires subcycle-separation != none and " *
            "subcycle-separation-engine == callback",
        ))
    end
    return nothing
end

"""Static data for Path-CBRP max-flow subtour separation (root loop or callback)."""
struct PathSubtourSepContext
    model::Model
    x
    y
    A::Arcs
    y_meta::Vector{Tuple{Int,Int}}
    out_idx::Dict{Int,Vector{Int}}
    depot::Int
    V::Vi
    y_at::Dict{Tuple{Int,Int},Int}
    blocks_at::Dict{Int,Vector{Int}}
    Vₘ::Dict{Int,Int}
    Vₘʳ::Dict{Int,Int}
    n::Int
    epsilon::Float64
    sep_mode::String
    M::Int
end

"""Map cluster index `b` and node `i` to flattened `y` index."""
function path_y_index_map(y_meta::Vector{Tuple{Int,Int}})::Dict{Tuple{Int,Int},Int}
    return Dict{Tuple{Int,Int},Int}((y_meta[k][1], y_meta[k][2]) => k for k in 1:length(y_meta))
end

"""`B(i)` = cluster indices `b` with `i` in block `b`."""
function path_blocks_at(
    y_meta::Vector{Tuple{Int,Int}},
)::Dict{Int,Vector{Int}}
    blocks_at::Dict{Int,Vector{Int}} = Dict{Int,Vector{Int}}()
    for k::Int in 1:length(y_meta)
        b::Int, i::Int = y_meta[k]
        push!(get!(Vector{Int}, blocks_at, i), b)
    end
    for i::Int in keys(blocks_at)
        blocks_at[i] = sort!(unique(blocks_at[i]))
    end
    return blocks_at
end

"""Arc indices `k` with `A[k] ∈ δ⁺(S)`."""
function path_out_arc_indices_S(
    S::Set{Int},
    A::Arcs,
    out_idx::Dict{Int,Vector{Int}},
)::Vector{Int}
    ks::Vector{Int} = Int[]
    for i::Int in S
        for k::Int in get(out_idx, i, Int[])
            last(A[k]) in S && continue
            push!(ks, k)
        end
    end
    return ks
end

"""
Build separation context for Path-CBRP subtour cuts.
"""
function buildPathSubtourSepContext(
    data::SBRPData,
    model::Model,
    app::Dict{String,Any},
    A::Arcs,
    y_meta::Vector{Tuple{Int,Int}},
    out_idx::Dict{Int,Vector{Int}},
    depot::Int,
)::PathSubtourSepContext
    V::Vi = Vi(collect(keys(data.D.V)))
    y_at::Dict{Tuple{Int,Int},Int} = path_y_index_map(y_meta)
    blocks_at::Dict{Int,Vector{Int}} = path_blocks_at(y_meta)
    Vₘ = Dict{Int,Int}(map((idx, i)::Tuple{Int,Int} -> i => idx, enumerate(V)))
    Vₘʳ = Dict{Int,Int}(map((idx, i)::Tuple{Int,Int} -> idx => i, enumerate(V)))
    return PathSubtourSepContext(
        model,
        model[:x],
        model[:y],
        A,
        y_meta,
        out_idx,
        depot,
        V,
        y_at,
        blocks_at,
        Vₘ,
        Vₘʳ,
        length(Vₘ),
        secMinViolation(app),
        get(app, "subcycle-separation", "all"),
        100_000,
    )
end

"""
Re-enable CPLEX presolve and aggregator for the final MIP after cut separation
(`getPathSubtourCuts` sets both to `0` during the LP separation loop).
"""
function restorePathCbrpMipPreprocessor!(model::Model)::Nothing
    set_optimizer_attribute(model, "CPXPARAM_LPMethod", 0)
    set_optimizer_attribute(model, "CPXPARAM_Preprocessing_Presolve", 1)
    set_optimizer_attribute(model, "CPXPARAM_Preprocessing_Aggregator", -1)
    return nothing
end

"""
Add Path SEC: `sum_{a in δ⁺(S)} x_a >= y_{ky_i} + y_{ky_j} - 1`.
"""
function addPathSubtourCut!(
    model::Model,
    A::Arcs,
    out_idx::Dict{Int,Vector{Int}},
    cut::PathSubtourCut,
)::Nothing
    S::Set{Int}, i::Int, j::Int, ky_i::Int, ky_j::Int = cut
    x = model[:x]
    y = model[:y]
    out_k::Vector{Int} = path_out_arc_indices_S(S, A, out_idx)
    c_ref = @constraint(
        model,
        sum(x[k] for k in out_k; init=0.0) >= y[ky_i] + y[ky_j] - 1,
    )
    set_name(c_ref, "c_path_sec[i=$(i),j=$(j),kyi=$(ky_i),kyj=$(ky_j),|S|=$(length(S))]")
    return nothing
end

"""
Submit Path SEC as a CPLEX user cut from a callback.
"""
function submitPathSubtourUserCut!(
    cb_data::CPLEX.CallbackContext,
    ctx::PathSubtourSepContext,
    cut::PathSubtourCut,
)::Nothing
    S::Set{Int}, _, _, ky_i::Int, ky_j::Int = cut
    out_k::Vector{Int} = path_out_arc_indices_S(S, ctx.A, ctx.out_idx)
    MOI.submit(
        ctx.model,
        MOI.UserCut(cb_data),
        @build_constraint(
            sum(ctx.x[k] for k in out_k; init=0.0) >= ctx.y[ky_i] + ctx.y[ky_j] - 1,
        ),
    )
    return nothing
end

"""
Submit a depot-rooted Path SEC `sum_{a in δ⁺(S)} x_a >= y_k` as a CPLEX lazy constraint
(integer candidates).
"""
function submitPathDepotSubtourLazyCut!(
    cb_data::CPLEX.CallbackContext,
    ctx::PathSubtourSepContext,
    cut::PathDepotSubtourCut,
)::Nothing
    S::Set{Int}, ky::Int = cut
    out_k::Vector{Int} = path_out_arc_indices_S(S, ctx.A, ctx.out_idx)
    MOI.submit(
        ctx.model,
        MOI.LazyConstraint(cb_data),
        @build_constraint(sum(ctx.x[k] for k in out_k; init=0.0) >= ctx.y[ky]),
    )
    return nothing
end

"""
\\brief Connected components of the integer `x`-support (arcs with `x > 0.5`), via iterative DFS.

Flow balance makes each weakly connected component Eulerian (hence strongly connected), so an
undirected search suffices.

\\param A Arc list indexed like `x_val`
\\param x_val Integer (within tolerance) arc values
\\return Map node → component id (only nodes incident to a support arc)
"""
function pathSupportComponents(A::Arcs, x_val::Dict{Int,Float64})::Dict{Int,Int}
    adj::Dict{Int,Vector{Int}} = Dict{Int,Vector{Int}}()
    for (k::Int, a::Arc) in enumerate(A)
        get(x_val, k, 0.0) > 0.5 || continue
        push!(get!(Vector{Int}, adj, first(a)), last(a))
        push!(get!(Vector{Int}, adj, last(a)), first(a))
    end
    comp::Dict{Int,Int} = Dict{Int,Int}()
    n_comp::Int = 0
    for s::Int in keys(adj)
        haskey(comp, s) && continue
        n_comp += 1
        comp[s] = n_comp
        stack::Vector{Int} = [s]
        while !isempty(stack)
            u::Int = pop!(stack)
            for v::Int in adj[u]
                haskey(comp, v) && continue
                comp[v] = n_comp
                push!(stack, v)
            end
        end
    end
    return comp
end

"""
\\brief Exact SEC separation at an integer candidate: one depot-rooted cut
`sum_{a in δ⁺(C)} x_a >= y_k` per serviced `y_k` whose node lies in a support component `C`
without the depot.

\\param ctx Path SEC separation context (uses `A`, `depot`, `y_meta`)
\\param x_val Integer arc values
\\param y_val Integer service values
\\return Violated cuts `(C, k)`; empty iff every serviced node is connected to the depot
"""
function findDisconnectedPathSubtourCuts(
    ctx::PathSubtourSepContext,
    x_val::Dict{Int,Float64},
    y_val::Dict{Int,Float64},
)::Vector{PathDepotSubtourCut}
    comp::Dict{Int,Int} = pathSupportComponents(ctx.A, x_val)
    depot_c::Int = get(comp, ctx.depot, 0)
    members::Dict{Int,Set{Int}} = Dict{Int,Set{Int}}()
    for (v::Int, c::Int) in comp
        push!(get!(Set{Int}, members, c), v)
    end
    cuts::Vector{PathDepotSubtourCut} = PathDepotSubtourCut[]
    for k::Int in 1:length(ctx.y_meta)
        get(y_val, k, 0.0) > 0.5 || continue
        c::Int = get(comp, ctx.y_meta[k][2], 0)
        (c == 0 || c == depot_c) && continue
        push!(cuts, (members[c], k))
    end
    return cuts
end

"""Largest `|v - round(v)|` over `x_val` and `y_val` (candidate integrality diagnostic)."""
function pathCandidateMaxFractionality(
    x_val::Dict{Int,Float64},
    y_val::Dict{Int,Float64},
)::Float64
    m::Float64 = 0.0
    for v in Iterators.flatten((values(x_val), values(y_val)))
        m = max(m, abs(v - round(v)))
    end
    return m
end

"""
Max-flow separation of violated Path-CBRP subtour cuts at a fractional `(x_val, y_val)`.
Max-flow (on truncated capacities) only proposes the sets `S`; a cut is reported iff its
violation on the exact `x_val`, `y_val` exceeds `ctx.epsilon`. Does not modify the model.
"""
function findViolatedPathSubtourCuts(
    ctx::PathSubtourSepContext,
    x_val::Dict{Int,Float64},
    y_val::Dict{Int,Float64},
)::Set{PathSubtourCut}
    new_cuts::Set{PathSubtourCut} = Set{PathSubtourCut}()
    max_violation::Float64 = 0.0

    g = SparseMaxFlowMinCut.ArcFlow[]
    A′_k::Vector{Int} = [k for k in 1:length(ctx.A) if x_val[k] > EPS]
    V′::Vi = Vi([
        i for i::Int in ctx.V if any(
            begin
                ky::Int = get(ctx.y_at, (b, i), 0)
                ky > 0 && get(y_val, ky, 0.0) > EPS
            end for b in get(ctx.blocks_at, i, Int[])
        )
    ])
    sources′::Vi = unique(vcat([ctx.depot], V′))

    for k::Int in A′_k
        a::Arc = ctx.A[k]
        push!(
            g,
            SparseMaxFlowMinCut.ArcFlow(
                ctx.Vₘ[first(a)],
                ctx.Vₘ[last(a)],
                trunc(floor(x_val[k], digits=5) * ctx.M),
            ),
        )
    end

    for source::Int in sources′
        for target::Int in V′
            source == target && continue

            _, _, set = SparseMaxFlowMinCut.find_maxflow_mincut(
                SparseMaxFlowMinCut.Graph(ctx.n, g),
                ctx.Vₘ[source],
                ctx.Vₘ[target],
            )

            set[ctx.Vₘ[target]] == 1 && continue

            S::Set{Int} = Set{Int}(
                map(i::Int -> ctx.Vₘʳ[i], filter(i::Int -> set[i] == 1, 1:ctx.n)),
            )
            # Capacities are truncated to 5 digits, so the max-flow value can undershoot
            # x(δ⁺(S)) by ~1e-5 per arc; violation is judged on the exact values.
            flow::Float64 = sum(
                get(x_val, k, 0.0) for k in path_out_arc_indices_S(S, ctx.A, ctx.out_idx);
                init=0.0,
            )

            for i::Int in S
                for j::Int in V′
                    j in S && continue
                    for b::Int in get(ctx.blocks_at, i, Int[])
                        ky_i::Int = get(ctx.y_at, (b, i), 0)
                        ky_i == 0 && continue
                        for b′::Int in get(ctx.blocks_at, j, Int[])
                            ky_j::Int = get(ctx.y_at, (b′, j), 0)
                            ky_j == 0 && continue
                            rhs::Float64 = y_val[ky_i] + y_val[ky_j] - 1.0
                            flow + ctx.epsilon >= rhs && continue

                            violation::Float64 = rhs - (flow + ctx.epsilon)
                            cut::PathSubtourCut = (S, i, j, ky_i, ky_j)

                            if ctx.sep_mode == "best"
                                if max_violation < violation
                                    empty!(new_cuts)
                                else
                                    continue
                                end
                            end

                            max_violation = max(violation, max_violation)
                            push!(new_cuts, cut)

                            if ctx.sep_mode == "first"
                                return new_cuts
                            end
                        end
                    end
                end
            end
        end
    end

    return new_cuts
end

"""
\\brief Canonical key of a Path SEC user cut; interns `S` in `stats.s_ids`.

The constraint depends only on `δ⁺(S)` and the unordered pair `{ky_i, ky_j}`, so `i`, `j` are dropped.

\\param stats Callback statistics holding the `S` intern table
\\param cut Path SEC `(S, i, j, ky_i, ky_j)`
\\return `(id(S), min(ky_i, ky_j), max(ky_i, ky_j))`
"""
function pathUserCutKey!(stats::PathSubtourCallbackStats, cut::PathSubtourCut)::PathUserCutKey
    S::Set{Int}, _, _, ky_i::Int, ky_j::Int = cut
    s_id::Int = get!(stats.s_ids, sort!(collect(S)), length(stats.s_ids) + 1)
    return (s_id, min(ky_i, ky_j), max(ky_i, ky_j))
end

"""
\\brief Violation `y_{ky_i} + y_{ky_j} - 1 - sum_{a in δ⁺(S)} x_a` of a Path SEC on untruncated values.

\\param ctx Path SEC separation context (uses `A`, `out_idx`)
\\param cut Path SEC `(S, i, j, ky_i, ky_j)`
\\param x_val Arc values
\\param y_val Service values
\\return Violation (positive iff the point violates the cut)
"""
function pathSubtourCutViolation(
    ctx::PathSubtourSepContext,
    cut::PathSubtourCut,
    x_val::Dict{Int,Float64},
    y_val::Dict{Int,Float64},
)::Float64
    S::Set{Int}, _, _, ky_i::Int, ky_j::Int = cut
    lhs::Float64 = sum(
        get(x_val, k, 0.0) for k in path_out_arc_indices_S(S, ctx.A, ctx.out_idx);
        init=0.0,
    )
    return get(y_val, ky_i, 0.0) + get(y_val, ky_j, 0.0) - 1.0 - lhs
end

"""
\\brief Filter one relaxation round of max-flow cuts against previously submitted keys and update stats.

Every found cut is checked against `stats.seen_user_cuts`; repeats are counted in `n_user_dup` and
dropped only when `stats.dedup_user_cuts`. Kept cuts are recorded, and their exact violation feeds the
`user_viol_*` stats (`<= ctx.epsilon` counts as non-violated). Does not touch `n_user_cuts`.

\\param stats Callback statistics (mutated)
\\param ctx Path SEC separation context
\\param cuts Cuts returned by `findViolatedPathSubtourCuts`
\\param x_val Arc values of the separated point
\\param y_val Service values of the separated point
\\return `(to_submit, round)` where `round` is a NamedTuple with `found`, `dup`, `submitted`,
        `nonviolated`, `viol_min`, `viol_mean`, `viol_max` for this round
"""
function selectPathUserCuts!(
    stats::PathSubtourCallbackStats,
    ctx::PathSubtourSepContext,
    cuts,
    x_val::Dict{Int,Float64},
    y_val::Dict{Int,Float64},
)
    to_submit::Vector{PathSubtourCut} = PathSubtourCut[]
    n_dup::Int = 0
    n_nonviolated::Int = 0
    viol_min::Float64 = Inf
    viol_max::Float64 = -Inf
    viol_sum::Float64 = 0.0
    for cut::PathSubtourCut in cuts
        key::PathUserCutKey = pathUserCutKey!(stats, cut)
        if key in stats.seen_user_cuts
            n_dup += 1
            stats.dedup_user_cuts && continue
        else
            push!(stats.seen_user_cuts, key)
        end
        v::Float64 = pathSubtourCutViolation(ctx, cut, x_val, y_val)
        v <= ctx.epsilon && (n_nonviolated += 1)
        viol_min = min(viol_min, v)
        viol_max = max(viol_max, v)
        viol_sum += v
        push!(to_submit, cut)
    end
    n_found::Int = length(cuts)
    stats.n_user_found += n_found
    stats.n_user_dup += n_dup
    stats.n_user_nonviolated += n_nonviolated
    n_found > 0 && (stats.n_user_rounds += 1)
    stats.user_viol_min = min(stats.user_viol_min, viol_min)
    stats.user_viol_max = max(stats.user_viol_max, viol_max)
    stats.user_viol_sum += viol_sum
    n_sub::Int = length(to_submit)
    return to_submit, (
        found=n_found,
        dup=n_dup,
        submitted=n_sub,
        nonviolated=n_nonviolated,
        viol_min=viol_min,
        viol_mean=n_sub > 0 ? viol_sum / n_sub : NaN,
        viol_max=viol_max,
    )
end

"""Mean exact violation of all submitted user cuts (`NaN` when none were submitted)."""
pathUserCutMeanViolation(stats::PathSubtourCallbackStats)::Float64 =
    stats.n_user_cuts > 0 ? stats.user_viol_sum / stats.n_user_cuts : NaN

"""Short scientific rendering for violation diagnostics."""
_viol_str(v::Float64)::String = isfinite(v) ? string(round(v; sigdigits=3)) : "NA"

"""Load `(x_val, y_val)` from a CPLEX generic callback after `load_callback_variable_primal`."""
function pathCallbackPrimalValues(
    cb_data::CPLEX.CallbackContext,
    ctx::PathSubtourSepContext,
)::Tuple{Dict{Int,Float64},Dict{Int,Float64}}
    x_val::Dict{Int,Float64} = Dict{Int,Float64}(
        k => callback_value(cb_data, ctx.x[k]) for k in 1:length(ctx.A)
    )
    y_val::Dict{Int,Float64} = Dict{Int,Float64}(
        k => callback_value(cb_data, ctx.y[k]) for k in 1:length(ctx.y_meta)
    )
    return x_val, y_val
end

"""
Log Path SEC callback submissions to stdout (interleaves with CPLEX `User` / `UserPurge` lines).

Disable with `PATH_CBRP_SEC_CALLBACK_LOG=0`.
"""
function pathSecCallbackLog!(
    context::String,
    n_added::Int,
    stats::PathSubtourCallbackStats,
)::Nothing
    get(ENV, "PATH_CBRP_SEC_CALLBACK_LOG", "1") == "0" && return nothing
    n_added == 0 && return nothing
    println(
        "[PathSEC] $(context): +$(n_added) cuts " *
        "(cum. user=$(stats.n_user_cuts), lazy=$(stats.n_lazy_cuts))",
    )
    flush(stdout)
    return nothing
end

"""
Log one RELAXATION round (found / duplicate / submitted cuts and exact violations).

Silent when separation found nothing; disable with `PATH_CBRP_SEC_CALLBACK_LOG=0`.
"""
function pathSecRelaxationLog!(rnd, stats::PathSubtourCallbackStats)::Nothing
    get(ENV, "PATH_CBRP_SEC_CALLBACK_LOG", "1") == "0" && return nothing
    rnd.found == 0 && return nothing
    println(
        "[PathSEC] RELAXATION user: +$(rnd.submitted) cuts " *
        "[found=$(rnd.found), dup=$(rnd.dup), nonViolated=$(rnd.nonviolated), " *
        "viol min/mean/max=$(_viol_str(rnd.viol_min))/$(_viol_str(rnd.viol_mean))/" *
        "$(_viol_str(rnd.viol_max))] " *
        "(cum. user=$(stats.n_user_cuts), found=$(stats.n_user_found), " *
        "dup=$(stats.n_user_dup), lazy=$(stats.n_lazy_cuts))",
    )
    flush(stdout)
    return nothing
end

"""True when all `x` and `y` callback values are (near) binary integers."""
function pathCallbackSolutionIsInteger(
    x_val::Dict{Int,Float64},
    y_val::Dict{Int,Float64},
)::Bool
    for v in values(x_val)
        abs(v - round(v)) > 1e-5 && return false
    end
    for v in values(y_val)
        abs(v - round(v)) > 1e-5 && return false
    end
    return true
end

"""
Register a single CPLEX callback: user cuts at LP relaxations (nodes up to
`stats.user_cut_max_depth`), lazy cuts at integer candidates (every node).

CPLEX.jl convention (see `CPLEX/test/MathOptInterface/MOI_callbacks.jl`):
- `CPX_CALLBACKCONTEXT_RELAXATION` → `MOI.UserCut` (max-flow `findViolatedPathSubtourCuts`)
- `CPX_CALLBACKCONTEXT_CANDIDATE` → `MOI.LazyConstraint` (integer points only; exact DFS
  `findDisconnectedPathSubtourCuts`, depot-rooted cuts)

Do not use `callback_node_status` here: for generic callbacks CPLEX.jl maps every
`CANDIDATE` to `CALLBACK_NODE_STATUS_INTEGER`, so fractional LP points were never cut.
"""
function registerPathSubtourSeparationCallback!(
    model::Model,
    ctx::PathSubtourSepContext,
    stats::PathSubtourCallbackStats,
)::Nothing
    function path_subtour_sep_cb(cb_data::CPLEX.CallbackContext, context_id::Clong)
        n_added::Int = 0
        sep_elapsed::Float64 = 0.0
        if context_id == CPLEX.CPX_CALLBACKCONTEXT_RELAXATION
            if !secUserCutDepthAllowed(stats.user_cut_max_depth, callbackNodeDepth(cb_data))
                stats.n_user_skipped_depth += 1
                return nothing
            end
            CPLEX.load_callback_variable_primal(cb_data, context_id)
            x_val::Dict{Int,Float64}, y_val::Dict{Int,Float64} =
                pathCallbackPrimalValues(cb_data, ctx)
            sep_round = nothing
            sep_elapsed = @elapsed begin
                cuts::Set{PathSubtourCut} =
                    findViolatedPathSubtourCuts(ctx, x_val, y_val)
                to_submit, sep_round = selectPathUserCuts!(stats, ctx, cuts, x_val, y_val)
                for cut in to_submit
                    submitPathSubtourUserCut!(cb_data, ctx, cut)
                    stats.n_user_cuts += 1
                end
            end
            stats.sep_time += sep_elapsed
            pathSecRelaxationLog!(sep_round, stats)
        elseif context_id == CPLEX.CPX_CALLBACKCONTEXT_CANDIDATE
            ispoint_p = Ref{CPLEX.CPXINT}()
            if CPLEX.CPXcallbackcandidateispoint(cb_data, ispoint_p) != 0 ||
               ispoint_p[] == 0
                return
            end
            CPLEX.load_callback_variable_primal(cb_data, context_id)
            x_val, y_val = pathCallbackPrimalValues(cb_data, ctx)
            pathCallbackSolutionIsInteger(x_val, y_val) || return nothing
            n_added = 0
            lazy_cuts::Vector{PathDepotSubtourCut} = PathDepotSubtourCut[]
            sep_elapsed = @elapsed begin
                lazy_cuts = findDisconnectedPathSubtourCuts(ctx, x_val, y_val)
                for cut in lazy_cuts
                    submitPathDepotSubtourLazyCut!(cb_data, ctx, cut)
                    stats.n_lazy_cuts += 1
                    n_added += 1
                end
            end
            stats.sep_time += sep_elapsed
            if !isempty(lazy_cuts)
                comp_sizes::Vector{Int} = sort!(unique(map(c -> length(first(c)), lazy_cuts)))
                pathSecCallbackLog!(
                    "CANDIDATE lazy [disconnected |C|=$(comp_sizes), " *
                    "maxFrac=$(pathCandidateMaxFractionality(x_val, y_val))]",
                    n_added,
                    stats,
                )
            end
        end
        return nothing
    end
    ctx_mask::UInt16 =
        CPLEX.CPX_CALLBACKCONTEXT_CANDIDATE | CPLEX.CPX_CALLBACKCONTEXT_RELAXATION
    MOI.set(model, CPLEX.CallbackFunction(ctx_mask), path_subtour_sep_cb)
    return nothing
end

"""Deprecated: use `registerPathSubtourSeparationCallback!`."""
registerPathSubtourUserCutCallback!(model, ctx, stats) =
    registerPathSubtourSeparationCallback!(model, ctx, stats)

"""
Max-flow separation of Path-CBRP subtour cuts (coupling `x` and `y_{b,i}`) via pre-MIP LP loop.
"""
function getPathSubtourCuts(
    data::SBRPData,
    model::Model,
    app::Dict{String,Any},
    A::Arcs,
    y_meta::Vector{Tuple{Int,Int}},
    out_idx::Dict{Int,Vector{Int}},
    depot::Int,
)::Set{PathSubtourCut}

    @debug "Path-CBRP: obtaining subtour cuts (y-coupled SEC)"

    ctx::PathSubtourSepContext =
        buildPathSubtourSepContext(data, model, app, A, y_meta, out_idx, depot)
    components::Set{PathSubtourCut} = Set{PathSubtourCut}()
    iteration::Int = 1

    set_optimizer_attribute(model, "CPXPARAM_LPMethod", 2)
    set_optimizer_attribute(model, "CPXPARAM_Advance", 1)

    while true
        _set_cplex_threads!(model, 1)
        @elapsed optimize!(model)

        if !in(termination_status(model), [MOI.OPTIMAL, MOI.TIME_LIMIT, MOI.ALMOST_INFEASIBLE])
            throw(InvalidStateException("Path-CBRP subtour separation: model could not be solved"))
        end

        y_val::Dict{Int,Float64} =
            Dict{Int,Float64}(k => value(ctx.y[k]) for k in 1:length(y_meta))
        x_val::Dict{Int,Float64} =
            Dict{Int,Float64}(k => value(ctx.x[k]) for k in 1:length(A))

        new_cuts::Set{PathSubtourCut} = findViolatedPathSubtourCuts(ctx, x_val, y_val)

        isempty(new_cuts) && break
        for cut in new_cuts
            addPathSubtourCut!(model, A, out_idx, cut)
        end
        union!(components, new_cuts)

        if iteration == 1
            set_optimizer_attribute(model, "CPXPARAM_Preprocessing_Presolve", 0)
            set_optimizer_attribute(model, "CPXPARAM_Preprocessing_Aggregator", 0)
        end
        iteration += 1
    end

    return components
end
