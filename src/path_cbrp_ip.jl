using Logging
using Printf

"""
Arc-indexed Path-CBRP MILP on a sparse Carlos digraph: compact arc MTZ (`w` per arc) plus global time budget.
Mirrors the C++ `PathCBRPIP` formulation in the parent COPAlgorithms repository (`src/directed/path_cbrp_model.cpp`).
"""
function runPathCbrpMipModel(data::SBRPData, app::Dict{String,Any})::Tuple{SBRPSolution,Dict{String,String}}

    depot::Int = depot_node(data)
    n_clusters::Int = num_clusters(data)
    A::Arcs = Arcs(sort(collect(data.D.A); by=a -> (first(a), last(a))))
    na::Int = length(A)
    if na == 0
        throw(ArgumentError("runPathCbrpMipModel: digraph has no arcs"))
    end

    out_idx::Dict{Int,Vector{Int}} = Dict{Int,Vector{Int}}()
    in_idx::Dict{Int,Vector{Int}} = Dict{Int,Vector{Int}}()
    for (k, a) in enumerate(A)
        push!(get!(Vector{Int}, out_idx, first(a)), k)
        push!(get!(Vector{Int}, in_idx, last(a)), k)
    end

    # Flatten y_{b,i} for every block vertex (Julia blocks are 1..n_clusters; no depot cluster in `B`).
    y_meta::Vector{Tuple{Int,Int}} = Tuple{Int,Int}[]
    for b::Int in 1:n_clusters
        for i::Int in get_cluster(data, b)
            push!(y_meta, (b, i))
        end
    end
    ny::Int = length(y_meta)

    info::Dict{String,String} = Dict{String,String}(
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
        "warmStartUsed" => "false",
        "warmStartFixed" => "false",
    )

    total_budget_sec::Float64 = _ip_base_cplex_seconds(parse(Float64, app["time-limit"]))
    phase1_wall_start::Float64 = Base.time()

    model::Model = direct_model(CPLEX.Optimizer())
    _set_cplex_threads!(model)
    set_time_limit_sec(model, max(1.0, total_budget_sec))

    #set_optimizer_attribute(model, "CPXPARAM_Preprocessing_Presolve", 0)

    @variable(model, x[1:na], Bin)
    @variable(model, y[1:ny], Bin)
    @variable(model, w[1:na], lower_bound = 0.0, upper_bound = data.T)

    @objective(model, Max, sum(k -> profit(data, y_meta[k][1]) * y[k], 1:ny))

    out_depot = get(out_idx, depot, Int[])
    in_depot = get(in_idx, depot, Int[])
    isempty(out_depot) && throw(ArgumentError("Path-CBRP: no outgoing arcs from depot"))
    isempty(in_depot) && throw(ArgumentError("Path-CBRP: no incoming arcs to depot"))
    @constraint(model, c_depot_leave_once, sum(x[k] for k in out_depot) == 1)
    @constraint(model, c_depot_enter_once, sum(x[k] for k in in_depot) == 1)

    verts_non_depot::Vector{Int} = sort!(collect(Int, filter(v::Int -> v != depot, keys(data.D.V))))
    @constraint(
        model,
        c_flow_balance[i=verts_non_depot],
        sum(x[k] for k in get(out_idx, i, Int[]); init=0.0) -
        sum(x[k] for k in get(in_idx, i, Int[]); init=0.0) == 0,
    )

    # At most one serviced vertex per block.
    block_indices::Dict{Int,Vector{Int}} = Dict{Int,Vector{Int}}()
    for k::Int in 1:ny
        b::Int = y_meta[k][1]
        push!(get!(Vector{Int}, block_indices, b), k)
    end
    block_ids_sorted::Vector{Int} = sort!(collect(Int, keys(block_indices)))
    @constraint(
        model,
        c_block_at_most_one[b=block_ids_sorted],
        sum(y[k] for k in block_indices[b]) <= 1,
    )


    # Outgoing tour arc if block served at node i.
    for k::Int in 1:ny
        _, i::Int = y_meta[k]
        o::Vector{Int} = get(out_idx, i, Int[])
        isempty(o) && throw(ArgumentError("Path-CBRP: node $i has no outgoing arcs but appears in a block"))
    end
    @constraint(
        model,
        c_y_requires_leave_arc[k=1:ny],
        sum(x[j] for j in get(out_idx, y_meta[k][2], Int[]); init=0.0) >= y[k],
    )

    #=
    @constraint(
        model,
        c_global_time_budget,
        sum(k -> arc_time(data, A[k]) * x[k], 1:na) +
        sum(k -> service_time(data, y_meta[k][1]) * y[k], 1:ny) <= data.T,
    )
        =#

    # Compact arc MTZ: for each non-depot i, a ∈ δ⁺(i), a' ∈ δ⁻(i):
    # w_{a'} >= w_a + x_a t_a - (2 - x_{a'} - x_a) T.
    Tlim::Float64 = data.T
    for i::Int in verts_non_depot
        for k_out::Int in get(out_idx, i, Int[])
            ta::Float64 = arc_time(data, A[k_out])
            for k_in::Int in get(in_idx, i, Int[])
                c_ref = @constraint(
                    model,
                    w[k_in] >=
                    w[k_out] + ta * x[k_out] - Tlim * (2 - x[k_in] - x[k_out]),
                )
                set_name(
                    c_ref,
                    "c_mtz_compact[i=$(i),ko=$(k_out),ki=$(k_in)]",
                )
            end
        end
    end
    @constraint(model, c_w_depot_leave_ub[k=out_depot], w[k] <= Tlim)

    fix_mode::Bool = pathCbrpFixWarmStartFlag(app)
    warm_xy = get(app, "path_cbrp_warm_xy", nothing)
    warm_sol = get(app, "warm_start_solution", nothing)
    if warm_xy !== nothing
        try
            x_on::Set{Arc}, y_on::Set{Tuple{Int,Int}} = warm_xy
            if fix_mode
                pathCbrpFixWarmFromXY!(model, A, y_meta, x_on, y_on)
                info["warmStartFixed"] = "true"
            else
                pathCbrpMipStartFromXY!(model, data, A, y_meta, x_on, y_on)
            end
            info["warmStartUsed"] = "true"
        catch e
            @warn "Path-CBRP MIP start skipped: $(sprint(showerror, e))"
            info["warmStartUsed"] = "false"
        end
    elseif warm_sol !== nothing
        try
            if fix_mode
                pathCbrpFixWarmFromSolution!(model, data, A, y_meta, warm_sol::SBRPSolution)
                info["warmStartFixed"] = "true"
            else
                pathCbrpMipStart!(model, data, A, y_meta, warm_sol::SBRPSolution)
            end
            info["warmStartUsed"] = "true"
        catch e
            @warn "Path-CBRP MIP start skipped: $(sprint(showerror, e))"
            info["warmStartUsed"] = "false"
        end
    end

    logPathCbrpWarmStartTimeAccounting!(data, app)

    info["phase1Time"] = string(Base.time() - phase1_wall_start)

    if get(app, "subcycle-separation", "none") != "none"
        unsetBinary(values(x))
        unsetBinary(values(y))
        sep_start::Float64 = Base.time()
        path_subtour_cuts::Set{PathSubtourCut} =
            getPathSubtourCuts(data, model, app, A, y_meta, out_idx, depot)
        restorePathCbrpMipPreprocessor!(model)
        info["maxFlowCutsTime"] = string(Base.time() - sep_start)
        info["maxFlowCuts"] = string(length(path_subtour_cuts))
        setBinary(values(x))
        setBinary(values(y))
    end

    info["solverTime"] = string(@elapsed optimize!(model))

    if !has_values(model)
        info["cost"] = "0.00"
        info["bestBound"] = _mip_best_bound_str(model)
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
    info["bestBound"] = _mip_best_bound_str(model)
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
        arc_idx::Union{Int,Nothing} = findfirst(a::Arc -> first(a) == cur, solution_arcs)
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
