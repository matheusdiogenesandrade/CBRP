using Logging
using Printf

"""
Arc-indexed Path-CBRP MILP on a sparse Carlos digraph: compact arc MTZ (`w` per arc) plus global time budget.
Mirrors the C++ `PathCBRPIP` formulation in the parent COPAlgorithms repository (`src/directed/path_cbrp_model.cpp`).
"""
function runPathCbrpMipModel(data::SBRPData, app::Dict{String, Any})::Tuple{SBRPSolution, Dict{String, String}}

    depot::Int = depot_node(data)
    n_clusters::Int = num_clusters(data)
    A::Arcs = Arcs(sort(collect(data.D.A); by = a -> (first(a), last(a))))
    na::Int = length(A)
    if na == 0
        throw(ArgumentError("runPathCbrpMipModel: digraph has no arcs"))
    end

    out_idx::Dict{Int, Vector{Int}} = Dict{Int, Vector{Int}}()
    in_idx::Dict{Int, Vector{Int}} = Dict{Int, Vector{Int}}()
    for (k, a) in enumerate(A)
        push!(get!(Vector{Int}, out_idx, first(a)), k)
        push!(get!(Vector{Int}, in_idx, last(a)), k)
    end

    # Flatten y_{b,i} for every block vertex (Julia blocks are 1..n_clusters; no depot cluster in `B`).
    y_meta::Vector{Tuple{Int, Int}} = Tuple{Int, Int}[]
    for b::Int in 1:n_clusters
        for i::Int in get_cluster(data, b)
            push!(y_meta, (b, i))
        end
    end
    ny::Int = length(y_meta)

    info::Dict{String, String} = Dict{String, String}(
        "intersectionCutsTime" => "0",
        "intersectionCuts1" => "0",
        "intersectionCuts2" => "0",
        "initialLPTime" => "0",
        "initialLP" => "N/A",
        "yLPTime" => "0",
        "yLP" => "N/A",
        "zLPTime" => "0",
        "zLP" => "N/A",
        "wLPTime" => "0",
        "wLP" => "N/A",
        "maxFlowLP" => "N/A",
        "maxFlowCuts" => "0",
        "maxFlowCutsTime" => "0",
        "phase1Time" => "0",
    )

    total_budget_sec::Float64 = _ip_base_cplex_seconds(parse(Float64, app["time-limit"]))
    phase1_wall_start::Float64 = Base.time()

    model::Model = direct_model(CPLEX.Optimizer())
    _set_cplex_threads!(model)
    set_time_limit_sec(model, max(1.0, total_budget_sec))

    @variable(model, x[1:na], Bin)
    @variable(model, y[1:ny], Bin)
    @variable(model, w[1:na], lower_bound = 0.0, upper_bound = data.T)

    @objective(model, Max, sum(k -> profit(data, y_meta[k][1]) * y[k], 1:ny))

    out_depot = get(out_idx, depot, Int[])
    in_depot = get(in_idx, depot, Int[])
    isempty(out_depot) && throw(ArgumentError("Path-CBRP: no outgoing arcs from depot"))
    isempty(in_depot) && throw(ArgumentError("Path-CBRP: no incoming arcs to depot"))
    @constraint(model, sum(x[k] for k in out_depot) == 1)
    @constraint(model, sum(x[k] for k in in_depot) == 1)

    V_all::Vi = Vi(collect(keys(data.D.V)))
    for i::Int in V_all
        if i == depot
            continue
        end
        o::Vector{Int} = get(out_idx, i, Int[])
        ine::Vector{Int} = get(in_idx, i, Int[])
        @constraint(model, sum(x[k] for k in o; init = 0.0) - sum(x[k] for k in ine; init = 0.0) == 0)
    end

    # At most one serviced vertex per block.
    block_indices::Dict{Int, Vector{Int}} = Dict{Int, Vector{Int}}()
    for k::Int in 1:ny
        b::Int = y_meta[k][1]
        push!(get!(Vector{Int}, block_indices, b), k)
    end
    for ks::Vector{Int} in values(block_indices)
        @constraint(model, sum(y[k] for k in ks) <= 1)
    end

    # Outgoing tour arc if block served at node i.
    for k::Int in 1:ny
        _, i::Int = y_meta[k]
        o::Vector{Int} = get(out_idx, i, Int[])
        isempty(o) && throw(ArgumentError("Path-CBRP: node $i has no outgoing arcs but appears in a block"))
        @constraint(model, sum(x[j] for j in o) >= y[k])
    end

    @constraint(
        model,
        sum(k -> arc_time(data, A[k]) * x[k], 1:na) +
        sum(k -> service_time(data, y_meta[k][1]) * y[k], 1:ny) <= data.T
    )

    # Compact arc MTZ: for each non-depot i, a ∈ δ⁺(i), a' ∈ δ⁻(i):
    # w_{a'} >= w_a + x_a t_a - (2 - x_{a'} - x_a) T.
    Tlim::Float64 = data.T
    for i::Int in V_all
        if i == depot
            continue
        end
        for k_out::Int in get(out_idx, i, Int[])
            ta::Float64 = arc_time(data, A[k_out])
            for k_in::Int in get(in_idx, i, Int[])
                @constraint(
                    model,
                    w[k_in] >=
                    w[k_out] + ta * x[k_out] - Tlim * (2 - x[k_in] - x[k_out]),
                )
            end
        end
    end
    for k::Int in out_depot
        @constraint(model, w[k] <= Tlim)
    end

    info["phase1Time"] = string(Base.time() - phase1_wall_start)
    info["solverTime"] = string(@elapsed optimize!(model))

    if !has_values(model)
        info["cost"] = "0.00"
        info["relativeGAP"] = "N/A"
        info["nodeCount"] = try
            string(node_count(model))
        catch
            "0"
        end
        info["noFeasibleSolution"] = "true"
        return SBRPSolution(Vi([depot]), VVi()), info
    end

    info["cost"] = @sprintf("%.2f", objective_value(model))
    info["relativeGAP"] = try
        string(relative_gap(model))
    catch
        "N/A"
    end
    info["nodeCount"] = try
        string(node_count(model))
    catch
        "0"
    end

    solution_arcs::Arcs = Arcs([A[k] for k in 1:na if value(x[k]) > 0.5])

    tour::Vi = Vi([depot])
    max_steps::Int = na + 5
    steps::Int = 0
    while steps < max_steps && !isempty(solution_arcs)
        steps += 1
        cur::Int = last(tour)
        arc_idx::Union{Int, Nothing} = findfirst(a::Arc -> first(a) == cur, solution_arcs)
        arc_idx === nothing && error("Path-CBRP: broken tour — arc not found leaving $(cur)")
        a::Arc = solution_arcs[arc_idx]
        push!(tour, last(a))
        deleteat!(solution_arcs, arc_idx)
        if last(tour) == depot && length(tour) > 1
            break
        end
    end
    if !isempty(solution_arcs)
        @warn "Path-CBRP: $(length(solution_arcs)) arc(s) remain after closing at depot (MTZ may still allow degeneracy); tour omits them."
    end

    solution_blocks::VVi = VVi()
    for k::Int in 1:ny
        if value(y[k]) > 0.5
            b::Int = y_meta[k][1]
            push!(solution_blocks, Vi(collect(Int, get_cluster(data, b))))
        end
    end

    return SBRPSolution(tour, solution_blocks), info
end
