using Logging
using CPLEX

"""Statistics collected by the complete-digraph SEC separation callback."""
mutable struct CompleteSubtourCallbackStats
    n_user_cuts::Int
    n_lazy_cuts::Int
    sep_time::Float64
end

"""
Require enabled SEC separation when using the callback engine on the complete model.
"""
function validateCompleteSubtourSeparation!(app::Dict{String,Any})::Nothing
    get(app, "subcycle-separation-engine", "root") != "callback" && return nothing
    get(app, "subcycle-separation", "none") == "none" &&
        throw(ArgumentError(
            "subcycle-separation-engine callback requires subcycle-separation != none",
        ))
    return nothing
end

"""Static data for complete-digraph max-flow subtour separation (root loop or callback)."""
struct CompleteSubtourSepContext
    model::Model
    x
    A::Arcs
    V::Vi
    Vₘ::Dict{Int,Int}
    Vₘʳ::Dict{Int,Int}
    n::Int
    epsilon::Float64
    sep_mode::String
    M::Int
end

function buildCompleteSubtourSepContext(
    data::SBRPData,
    model::Model,
    app::Dict{String,Any},
)::CompleteSubtourSepContext
    V::Vi = Vi(collect(keys(data.D.V)))
    Vₘ = Dict{Int,Int}(map((idx, i)::Tuple{Int,Int} -> i => idx, enumerate(V)))
    Vₘʳ = Dict{Int,Int}(map((idx, i)::Tuple{Int,Int} -> idx => i, enumerate(V)))
    return CompleteSubtourSepContext(
        model,
        model[:x],
        data.D.A,
        V,
        Vₘ,
        Vₘʳ,
        length(Vₘ),
        1e-2,
        get(app, "subcycle-separation", "all"),
        100_000,
    )
end

"""
Max-flow separation of violated complete-digraph subtour cuts at `(x_val, z_val, w_val)`.
Does not modify the model.
"""
function findViolatedCompleteSubtourCuts(
    ctx::CompleteSubtourSepContext,
    x_val::ArcCostMap,
    z_val::Dict{Int,Float64},
    w_val::Dict{Int,Float64},
)::Set{Tuple{Arcs,Arcs}}
    new_components::Set{Tuple{Arcs,Arcs}} = Set{Tuple{Arcs,Arcs}}()
    max_violation::Float64 = 0.0

    A′::Arcs = filter(a::Arc -> x_val[a] > EPS, ctx.A)
    V′::Vi = filter(i::Int -> get(z_val, i, 0.0) > EPS, ctx.V)
    depots′::Vi = filter(i::Int -> get(w_val, i, 0.0) > EPS, ctx.V)

    g = SparseMaxFlowMinCut.ArcFlow[]
    for a::Arc in A′
        i::Int, j::Int = first(a), last(a)
        push!(
            g,
            SparseMaxFlowMinCut.ArcFlow(
                ctx.Vₘ[i],
                ctx.Vₘ[j],
                trunc(floor(x_val[a], digits=5) * ctx.M),
            ),
        )
    end

    for source::Int in depots′
        for target::Int in V′
            source == target && continue

            maxFlow::Float64, flows, set = SparseMaxFlowMinCut.find_maxflow_mincut(
                SparseMaxFlowMinCut.Graph(ctx.n, g),
                ctx.Vₘ[source],
                ctx.Vₘ[target],
            )
            flow::Float64 = maxFlow / ctx.M

            set[ctx.Vₘ[target]] == 1 && continue
            flow + ctx.epsilon >= get(z_val, source, 0.0) + get(z_val, target, 0.0) - 1 &&
                continue

            S::Si = Si(
                map(
                    i::Int -> ctx.Vₘʳ[i],
                    filter(i::Int -> set[i] == 1, 1:ctx.n),
                ),
            )
            Aₛ::Arcs = δ⁺(ctx.A, S)
            Aᵢ::Arcs = union(δ⁺(A′, source), δ⁺(A′, target))
            violation::Float64 =
                get(z_val, source, 0.0) + get(z_val, target, 0.0) - 1 - (flow + ctx.epsilon)

            if ctx.sep_mode == "best"
                if max_violation < violation
                    empty!(new_components)
                else
                    continue
                end
            end

            max_violation = max(violation, max_violation)
            push!(new_components, (Aₛ, Aᵢ))

            if ctx.sep_mode == "first"
                return new_components
            end
        end

        if ctx.sep_mode == "first" && !isempty(new_components)
            return new_components
        end
    end

    return new_components
end

"""Submit complete subtour SEC as a CPLEX user cut from a callback."""
function submitCompleteSubtourUserCut!(
    cb_data::CPLEX.CallbackContext,
    ctx::CompleteSubtourSepContext,
    lhs_arcs::Arcs,
    rhs_arcs::Arcs,
)::Nothing
    MOI.submit(
        ctx.model,
        MOI.UserCut(cb_data),
        @build_constraint(
            sum(ctx.x[a] for a in lhs_arcs; init=0.0) >=
            sum(ctx.x[a] for a in rhs_arcs; init=0.0) - 1,
        ),
    )
    return nothing
end

"""Submit complete subtour SEC as a CPLEX lazy constraint from a callback."""
function submitCompleteSubtourLazyCut!(
    cb_data::CPLEX.CallbackContext,
    ctx::CompleteSubtourSepContext,
    lhs_arcs::Arcs,
    rhs_arcs::Arcs,
)::Nothing
    MOI.submit(
        ctx.model,
        MOI.LazyConstraint(cb_data),
        @build_constraint(
            sum(ctx.x[a] for a in lhs_arcs; init=0.0) >=
            sum(ctx.x[a] for a in rhs_arcs; init=0.0) - 1,
        ),
    )
    return nothing
end

"""Load `(x_val, z_val, w_val)` from a CPLEX generic callback after `load_callback_variable_primal`."""
function completeCallbackPrimalValues(
    cb_data::CPLEX.CallbackContext,
    ctx::CompleteSubtourSepContext,
)::Tuple{ArcCostMap,Dict{Int,Float64},Dict{Int,Float64}}
    x_val::ArcCostMap = ArcCostMap(
        a => callback_value(cb_data, ctx.x[a]) for a in ctx.A
    )
    z_val::Dict{Int,Float64} = Dict{Int,Float64}(
        i => callback_value(cb_data, ctx.model[:z][i]) for i in ctx.V
    )
    w_val::Dict{Int,Float64} = Dict{Int,Float64}(
        i => callback_value(cb_data, ctx.model[:w][i]) for i in ctx.V
    )
    return x_val, z_val, w_val
end

"""
Log complete SEC callback submissions to stdout (interleaves with CPLEX `User` / `UserPurge` lines).

Disable with `COMPLETE_SEC_CALLBACK_LOG=0`.
"""
function completeSecCallbackLog!(
    context::String,
    n_added::Int,
    stats::CompleteSubtourCallbackStats,
)::Nothing
    get(ENV, "COMPLETE_SEC_CALLBACK_LOG", "1") == "0" && return nothing
    n_added == 0 && return nothing
    println(
        "[CompleteSEC] $(context): +$(n_added) cuts " *
        "(cum. user=$(stats.n_user_cuts), lazy=$(stats.n_lazy_cuts))",
    )
    flush(stdout)
    return nothing
end

"""True when all `x`, `z`, and `w` callback values are (near) binary integers."""
function completeCallbackSolutionIsInteger(
    x_val::ArcCostMap,
    z_val::Dict{Int,Float64},
    w_val::Dict{Int,Float64},
)::Bool
    for v in values(x_val)
        abs(v - round(v)) > 1e-5 && return false
    end
    for v in values(z_val)
        abs(v - round(v)) > 1e-5 && return false
    end
    for v in values(w_val)
        abs(v - round(v)) > 1e-5 && return false
    end
    return true
end

"""
Register a single CPLEX callback: user cuts at LP relaxations, lazy cuts at integer candidates.
"""
function registerCompleteSubtourSeparationCallback!(
    model::Model,
    ctx::CompleteSubtourSepContext,
    stats::CompleteSubtourCallbackStats,
)::Nothing
    function complete_subtour_sep_cb(cb_data::CPLEX.CallbackContext, context_id::Clong)
        n_added::Int = 0
        sep_elapsed::Float64 = 0.0
        if context_id == CPLEX.CPX_CALLBACKCONTEXT_RELAXATION
            CPLEX.load_callback_variable_primal(cb_data, context_id)
            x_val::ArcCostMap, z_val::Dict{Int,Float64}, w_val::Dict{Int,Float64} =
                completeCallbackPrimalValues(cb_data, ctx)
            n_added = 0
            sep_elapsed = @elapsed begin
                cuts::Set{Tuple{Arcs,Arcs}} =
                    findViolatedCompleteSubtourCuts(ctx, x_val, z_val, w_val)
                for (lhs_arcs::Arcs, rhs_arcs::Arcs) in cuts
                    submitCompleteSubtourUserCut!(cb_data, ctx, lhs_arcs, rhs_arcs)
                    stats.n_user_cuts += 1
                    n_added += 1
                end
            end
            stats.sep_time += sep_elapsed
            completeSecCallbackLog!("RELAXATION user", n_added, stats)
        elseif context_id == CPLEX.CPX_CALLBACKCONTEXT_CANDIDATE
            ispoint_p = Ref{CPLEX.CPXINT}()
            if CPLEX.CPXcallbackcandidateispoint(cb_data, ispoint_p) != 0 ||
               ispoint_p[] == 0
                return
            end
            CPLEX.load_callback_variable_primal(cb_data, context_id)
            x_val, z_val, w_val = completeCallbackPrimalValues(cb_data, ctx)
            completeCallbackSolutionIsInteger(x_val, z_val, w_val) || return nothing
            n_added = 0
            sep_elapsed = @elapsed begin
                cuts = findViolatedCompleteSubtourCuts(ctx, x_val, z_val, w_val)
                for (lhs_arcs, rhs_arcs) in cuts
                    submitCompleteSubtourLazyCut!(cb_data, ctx, lhs_arcs, rhs_arcs)
                    stats.n_lazy_cuts += 1
                    n_added += 1
                end
            end
            stats.sep_time += sep_elapsed
            completeSecCallbackLog!("CANDIDATE lazy", n_added, stats)
        end
        return nothing
    end
    ctx_mask::UInt16 =
        CPLEX.CPX_CALLBACKCONTEXT_CANDIDATE | CPLEX.CPX_CALLBACKCONTEXT_RELAXATION
    MOI.set(model, CPLEX.CallbackFunction(ctx_mask), complete_subtour_sep_cb)
    return nothing
end
