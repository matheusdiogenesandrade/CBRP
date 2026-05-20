using Logging
using CPLEX

"""Path-CBRP SEC: (S, i, j, k_yi, k_yj) with y_meta[k]=(b,i)."""
const PathSubtourCut = Tuple{Set{Int},Int,Int,Int,Int}

"""Statistics collected by the Path SEC separation callback (user + lazy cuts)."""
mutable struct PathSubtourCallbackStats
    n_user_cuts::Int
    n_lazy_cuts::Int
    sep_time::Float64
end

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
        1e-6,
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
Submit Path SEC as a CPLEX lazy constraint from a callback (integer candidates).
"""
function submitPathSubtourLazyCut!(
    cb_data::CPLEX.CallbackContext,
    ctx::PathSubtourSepContext,
    cut::PathSubtourCut,
)::Nothing
    S::Set{Int}, _, _, ky_i::Int, ky_j::Int = cut
    out_k::Vector{Int} = path_out_arc_indices_S(S, ctx.A, ctx.out_idx)
    MOI.submit(
        ctx.model,
        MOI.LazyConstraint(cb_data),
        @build_constraint(
            sum(ctx.x[k] for k in out_k; init=0.0) >= ctx.y[ky_i] + ctx.y[ky_j] - 1,
        ),
    )
    return nothing
end

"""
Max-flow separation of violated Path-CBRP subtour cuts at a fractional `(x_val, y_val)`.
Does not modify the model.
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

            maxFlow::Float64, flows, set = SparseMaxFlowMinCut.find_maxflow_mincut(
                SparseMaxFlowMinCut.Graph(ctx.n, g),
                ctx.Vₘ[source],
                ctx.Vₘ[target],
            )
            flow::Float64 = maxFlow / ctx.M

            set[ctx.Vₘ[target]] == 1 && continue

            S::Set{Int} = Set{Int}(
                map(i::Int -> ctx.Vₘʳ[i], filter(i::Int -> set[i] == 1, 1:ctx.n)),
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
Register a single CPLEX callback: user cuts at LP relaxations, lazy cuts at integer candidates.

CPLEX.jl convention (see `CPLEX/test/MathOptInterface/MOI_callbacks.jl`):
- `CPX_CALLBACKCONTEXT_RELAXATION` → `MOI.UserCut`
- `CPX_CALLBACKCONTEXT_CANDIDATE` → `MOI.LazyConstraint` (integer points only)

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
            CPLEX.load_callback_variable_primal(cb_data, context_id)
            x_val::Dict{Int,Float64}, y_val::Dict{Int,Float64} =
                pathCallbackPrimalValues(cb_data, ctx)
            n_added = 0
            sep_elapsed = @elapsed begin
                cuts::Set{PathSubtourCut} =
                    findViolatedPathSubtourCuts(ctx, x_val, y_val)
                for cut in cuts
                    submitPathSubtourUserCut!(cb_data, ctx, cut)
                    stats.n_user_cuts += 1
                    n_added += 1
                end
            end
            stats.sep_time += sep_elapsed
            pathSecCallbackLog!("RELAXATION user", n_added, stats)
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
            sep_elapsed = @elapsed begin
                cuts = findViolatedPathSubtourCuts(ctx, x_val, y_val)
                for cut in cuts
                    submitPathSubtourLazyCut!(cb_data, ctx, cut)
                    stats.n_lazy_cuts += 1
                    n_added += 1
                end
            end
            stats.sep_time += sep_elapsed
            pathSecCallbackLog!("CANDIDATE lazy", n_added, stats)
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
