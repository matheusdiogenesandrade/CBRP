using Logging

"""Path-CBRP SEC: (S, i, j, k_yi, k_yj) with y_meta[k]=(b,i)."""
const PathSubtourCut = Tuple{Set{Int},Int,Int,Int,Int}

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
Max-flow separation of Path-CBRP subtour cuts (coupling `x` and `y_{b,i}`).
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

    x = model[:x]
    y = model[:y]
    V::Vi = collect(keys(data.D.V))
    y_at::Dict{Tuple{Int,Int},Int} = path_y_index_map(y_meta)
    blocks_at::Dict{Int,Vector{Int}} = path_blocks_at(y_meta)

    components::Set{PathSubtourCut} = Set{PathSubtourCut}()
    Vₘ = Dict{Int,Int}(map((idx, i)::Tuple{Int,Int} -> i => idx, enumerate(V)))
    Vₘʳ = Dict{Int,Int}(map((idx, i)::Tuple{Int,Int} -> idx => i, enumerate(V)))
    n::Int = length(Vₘ)
    iteration::Int = 1
    epsilon::Float64 = 1e-2
    sep_mode::String = get(app, "subcycle-separation", "all")

    set_optimizer_attribute(model, "CPXPARAM_LPMethod", 2)
    set_optimizer_attribute(model, "CPXPARAM_Advance", 1)

    while true
        _set_cplex_threads!(model, 1)
        @elapsed optimize!(model)

        if !in(termination_status(model), [MOI.OPTIMAL, MOI.TIME_LIMIT, MOI.ALMOST_INFEASIBLE])
            throw(InvalidStateException("Path-CBRP subtour separation: model could not be solved"))
        end

        y_val::Dict{Int,Float64} = Dict{Int,Float64}(k => value(y[k]) for k in 1:length(y_meta))
        x_val::Dict{Int,Float64} = Dict{Int,Float64}(k => value(x[k]) for k in 1:length(A))

        g = SparseMaxFlowMinCut.ArcFlow[]
        M::Int = 100_000
        new_cuts::Set{PathSubtourCut} = Set{PathSubtourCut}()
        max_violation::Float64 = 0.0

        A′_k::Vector{Int} = [k for k in 1:length(A) if x_val[k] > EPS]
        V′::Vi = Vi([
            i for i::Int in V if any(
                begin
                    ky::Int = get(y_at, (b, i), 0)
                    ky > 0 && get(y_val, ky, 0.0) > EPS
                end for b in get(blocks_at, i, Int[])
            )
        ])
        sources′::Vi = unique(vcat([depot], V′))

        for k::Int in A′_k
            a::Arc = A[k]
            push!(
                g,
                SparseMaxFlowMinCut.ArcFlow(
                    Vₘ[first(a)],
                    Vₘ[last(a)],
                    trunc(floor(x_val[k], digits=5) * M),
                ),
            )
        end

        for source::Int in sources′
            for target::Int in V′
                source == target && continue

                maxFlow::Float64, flows, set = SparseMaxFlowMinCut.find_maxflow_mincut(
                    SparseMaxFlowMinCut.Graph(n, g),
                    Vₘ[source],
                    Vₘ[target],
                )
                flow::Float64 = maxFlow / M

                set[Vₘ[target]] == 1 && continue

                S::Set{Int} = Set{Int}(map(i::Int -> Vₘʳ[i], filter(i::Int -> set[i] == 1, 1:n)))

                for i::Int in S
                    for j::Int in V′
                        j in S && continue
                        for b::Int in get(blocks_at, i, Int[])
                            ky_i::Int = get(y_at, (b, i), 0)
                            ky_i == 0 && continue
                            for b′::Int in get(blocks_at, j, Int[])
                                ky_j::Int = get(y_at, (b′, j), 0)
                                ky_j == 0 && continue
                                rhs::Float64 = y_val[ky_i] + y_val[ky_j] - 1.0
                                flow + epsilon >= rhs && continue

                                violation::Float64 = rhs - (flow + epsilon)
                                cut::PathSubtourCut = (S, i, j, ky_i, ky_j)

                                if sep_mode == "best"
                                    if max_violation < violation
                                        empty!(new_cuts)
                                    else
                                        continue
                                    end
                                end

                                max_violation = max(violation, max_violation)
                                push!(new_cuts, cut)
                                addPathSubtourCut!(model, A, out_idx, cut)

                                if sep_mode == "first"
                                    break
                                end
                            end
                            sep_mode == "first" && !isempty(new_cuts) && break
                        end
                        sep_mode == "first" && !isempty(new_cuts) && break
                    end
                    sep_mode == "first" && !isempty(new_cuts) && break
                end

                sep_mode == "first" && !isempty(new_cuts) && break
            end
            sep_mode == "first" && !isempty(new_cuts) && break
        end

        isempty(new_cuts) && break
        union!(components, new_cuts)

        if iteration == 1
            set_optimizer_attribute(model, "CPXPARAM_Preprocessing_Presolve", 0)
            set_optimizer_attribute(model, "CPXPARAM_Preprocessing_Aggregator", 0)
        end
        iteration += 1
    end

    return components
end
