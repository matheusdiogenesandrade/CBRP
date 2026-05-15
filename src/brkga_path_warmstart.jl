using JuMP
using Logging

"""
    _expand_block_node_chain(paths, chain)::Vi

Stitch compact block-node sequence `chain` into a vertex list using metric-closure
interior vertices stored in `paths` (same contract as `retrieveOriginalDigraphSolution`
in `run.jl`).
"""
function _expand_block_node_chain(paths::Dict{Arc,Vi}, chain::Vi)::Vi
    n = length(chain)
    n == 0 && return Vi()
    n == 1 && return Vi([chain[1]])
    out::Vi = Vi()
    for i::Int in 1:(n-1)
        u::Int = chain[i]
        v::Int = chain[i+1]
        a::Arc = Arc(u, v)
        haskey(paths, a) || throw(ArgumentError("Missing path interiors for compact arc $a"))
        append!(out, chain[i])
        append!(out, paths[a])
    end
    push!(out, chain[n])
    return out
end

"""
    _strip_depot_ends!(seq::Vi, depot::Int)::Nothing

Remove leading/trailing depot occurrences from `seq` (mutates `seq`).
"""
function _strip_depot_ends!(seq::Vi, depot::Int)::Nothing
    while !isempty(seq) && first(seq) == depot
        popfirst!(seq)
    end
    while !isempty(seq) && last(seq) == depot
        pop!(seq)
    end
    return nothing
end

"""
    compactToSparseSbrpSolution(data_sparse, paths, sol_compact)::SBRPSolution

Expand a BRKGA / compact-graph `SBRPSolution` into a closed walk on the **sparse**
Carlos digraph `data_sparse`, using shortest-path interiors from `paths` produced
by metric closure (`createCompleteDigraph`).

`paths` must be the dict returned together with the compact instance for the **same**
instance file and reader flags as `data_sparse`.

Throws `ArgumentError` if depot legs or a compact arc cannot be realized on the
sparse digraph, or if `checkSBRPSolution` fails.
"""
function compactToSparseSbrpSolution(
    data_sparse::SBRPData,
    paths::Dict{Arc,Vi},
    sol_compact::SBRPSolution,
)::SBRPSolution
    depot::Int = data_sparse.depot
    seq::Vi = Vi(sol_compact.tour)
    _strip_depot_ends!(seq, depot)

    if isempty(seq)
        sol::SBRPSolution = SBRPSolution(Vi([depot]), VVi(sol_compact.B))
        checkSBRPSolution(data_sparse, sol)
        return sol
    end

    mid::Vi = _expand_block_node_chain(paths, seq)

    full::Vi = Vi([depot])
    a0::Arc = Arc(depot, first(mid))
    haskey(data_sparse.D.distance, a0) ||
        throw(ArgumentError("Missing sparse arc $a0 for depot→first expanded node"))
    append!(full, mid)
    a1::Arc = Arc(last(mid), depot)
    haskey(data_sparse.D.distance, a1) ||
        throw(ArgumentError("Missing sparse arc $a1 for last expanded node→depot"))
    push!(full, depot)

    sol = SBRPSolution(full, VVi(sol_compact.B))
    checkSBRPSolution(data_sparse, sol)
    return sol
end

"""
    pathCbrpMipStart!(model, data, A, y_meta, sol)::Nothing

JuMP MIP-start hints for Path-CBRP (`x` arc tour, `y` block-vertex selection, `w` arc MTZ).
Does not alter constraints or objective. All binaries are explicitly initialized.
"""
function pathCbrpMipStart!(
    model::Model,
    data::SBRPData,
    A::Arcs,
    y_meta::Vector{Tuple{Int,Int}},
    sol::SBRPSolution,
)::Nothing
    na::Int = length(A)
    ny::Int = length(y_meta)
    arc_idx::Dict{Arc,Int} = Dict{Arc,Int}(a => k for (k, a) in enumerate(A))
    x = model[:x]
    y = model[:y]
    w = model[:w]
    tour::Vi = sol.tour

    for k::Int in 1:na
        set_start_value(x[k], 0.0)
        set_start_value(w[k], 0.0)
    end
    for k::Int in 1:ny
        set_start_value(y[k], 0.0)
    end

    for i::Int in 1:(length(tour)-1)
        a::Arc = Arc(tour[i], tour[i+1])
        haskey(arc_idx, a) || throw(ArgumentError("Warm-start tour arc $a not in sparse arc list"))
        set_start_value(x[arc_idx[a]], 1.0)
    end

    served_y::Set{Tuple{Int,Int}} = Set{Tuple{Int,Int}}()
    for b::Int in 1:num_clusters(data)
        blk::Vi = get_cluster(data, b)
        any(cand::Vi -> Set(cand) == Set(blk), sol.B) || continue
        v_pick::Union{Nothing,Int} = nothing
        for node::Int in tour
            if node in blk
                v_pick = node
                break
            end
        end
        v_pick === nothing &&
            throw(ArgumentError("Warm start: served block $b has no tour vertex in that block"))
        push!(served_y, (b, v_pick))
    end

    for k::Int in 1:ny
        if y_meta[k] in served_y
            set_start_value(y[k], 1.0)
        end
    end

    Tlim::Float64 = data.T
    cum::Float64 = 0.0
    for i::Int in 1:(length(tour)-1)
        a::Arc = Arc(tour[i], tour[i+1])
        kx::Int = arc_idx[a]
        set_start_value(w[kx], min(cum, Tlim))
        cum += arc_time(data, a)
        v::Int = tour[i+1]
        for (b::Int, vv::Int) in served_y
            vv == v && (cum += service_time(data, b))
        end
        cum = min(cum, Tlim)
    end

    return nothing
end

"""
    _tourFromWarmArcs(depot, x_on, max_steps)::Union{Vi,Nothing}

Trace a depot-closed walk using arcs in `x_on` (same idea as Path-CBRP solution extraction).
"""
function _tourFromWarmArcs(depot::Int, x_on::Set{Arc}, max_steps::Int)::Union{Vi,Nothing}
    tour::Vi = Vi([depot])
    remaining::Arcs = Arcs(collect(x_on))
    steps::Int = 0
    while steps < max_steps && !isempty(remaining)
        steps += 1
        cur::Int = last(tour)
        idx::Union{Int,Nothing} = findfirst(a::Arc -> first(a) == cur, remaining)
        idx === nothing && return nothing
        a::Arc = remaining[idx]
        push!(tour, last(a))
        deleteat!(remaining, idx)
        if last(tour) == depot && length(tour) > 1
            return tour
        end
    end
    return nothing
end

"""
    logPathCbrpWarmStartTimeAccounting!(data, app)::Nothing

Log travel and service minutes implied by the Path-CBRP warm payload (`warm_start_solution` or
`path_cbrp_warm_xy`), using the same `time` / `blockTime` formulas as `c_global_time_budget`.
"""
function logPathCbrpWarmStartTimeAccounting!(data::SBRPData, app::Dict{String,Any})::Nothing
    warm_sol = get(app, "warm_start_solution", nothing)
    warm_xy = get(app, "path_cbrp_warm_xy", nothing)
    if warm_sol === nothing && warm_xy === nothing
        return nothing
    end
    Tlim::Float64 = data.T
    tol::Float64 = 1e-5
    if warm_sol !== nothing
        sol::SBRPSolution = warm_sol::SBRPSolution
        travel = tourDistance(data, sol.tour) / NORMAL_SPEED

        service = sum(blk -> blockTime(data, blk), sol.B; init=0.0)

        total = travel + service
        chk = tourTime(data, sol)
        viol = total > Tlim + tol
        if abs(chk - total) > max(1e-4, tol * max(1.0, abs(chk)))
            @warn "Path-CBRP warm SBRPSolution: travel+service differs from tourTime()" travel_plus_service = total tourTime_fn = chk delta = abs(chk - total)
        end
        @info "Path-CBRP warm-start time (SBRPSolution)" travel_minutes = travel service_minutes = service total_minutes = total budget_minutes = Tlim slack_minutes = (Tlim - total) exceeds_budget = viol num_serviced_blocks = length(sol.B)
    else
        x_on::Set{Arc}, y_on::Set{Tuple{Int,Int}} = warm_xy
        travel = sum(a -> floor(time(data, a)), x_on; init=0.0)

        bs::Set{Int} = Set{Int}()
        for pr::Tuple{Int,Int} in y_on
            push!(bs, pr[1])
        end
        service = sum(b -> service_time(data, b), bs; init=0.0)

        total = travel + service
        viol = total > Tlim + tol
        depot::Int = depot_node(data)
        na_est::Int = length(data.D.A)
        tour_tr::Union{Vi,Nothing} = _tourFromWarmArcs(depot, x_on, na_est + 5)
        if tour_tr !== nothing
            traced = tourDistance(data, tour_tr) / NORMAL_SPEED
            if abs(traced - travel) > max(1e-3, 1e-4 * max(1.0, travel))
                @warn "Path-CBRP warm X/Y: sum(t_a) over X arcs differs from depot tour traced along X (multiset / disconnected X?)" sum_arc_times_from_X = travel traced_depot_tour_minutes = traced
            end
        else
            @warn "Path-CBRP warm X/Y: could not trace depot-closed walk from X (reported travel is still sum of t_a over arcs in X)"
        end
        @info "Path-CBRP warm-start time (path-cbrp-mtz X/Y)" travel_minutes = travel service_minutes = service total_minutes = total budget_minutes = Tlim slack_minutes = (Tlim - total) exceeds_budget = viol num_serviced_blocks = length(bs)
    end
    return nothing
end

"""
    pathCbrpMipStartFromXY!(model, data, A, y_meta, x_on, y_on)::Nothing

MIP-start hints from explicit sparse `x` arcs and `y` assignments `(b, i)` (Julia block index).
"""
function pathCbrpMipStartFromXY!(
    model::Model,
    data::SBRPData,
    A::Arcs,
    y_meta::Vector{Tuple{Int,Int}},
    x_on::Set{Arc},
    y_on::Set{Tuple{Int,Int}},
)::Nothing
    na::Int = length(A)
    ny::Int = length(y_meta)
    arc_idx::Dict{Arc,Int} = Dict{Arc,Int}(a => k for (k, a) in enumerate(A))
    x = model[:x]
    y = model[:y]
    w = model[:w]

    for k::Int in 1:na
        set_start_value(x[k], 0.0)
        set_start_value(w[k], 0.0)
    end
    for k::Int in 1:ny
        set_start_value(y[k], 0.0)
    end

    for a::Arc in x_on
        haskey(arc_idx, a) || throw(ArgumentError("Warm-start arc $a not in sparse arc list"))
        set_start_value(x[arc_idx[a]], 1.0)
    end

    y_meta_set::Set{Tuple{Int,Int}} = Set(y_meta)
    for pair::Tuple{Int,Int} in y_on
        pair in y_meta_set || throw(ArgumentError("Warm-start y assignment $pair not in y_meta"))
        idx::Union{Int,Nothing} = findfirst(k::Int -> y_meta[k] == pair, 1:ny)
        idx === nothing && throw(ArgumentError("Warm-start y assignment $pair not in y_meta"))
        set_start_value(y[idx], 1.0)
    end

    depot::Int = depot_node(data)
    tour::Union{Vi,Nothing} = _tourFromWarmArcs(depot, x_on, na + 5)
    if tour === nothing
        @warn "Path-CBRP warm start: could not trace tour from X arcs; skipping w hints"
        return nothing
    end

    Tlim::Float64 = data.T
    cum::Float64 = 0.0
    for i::Int in 1:(length(tour)-1)
        a::Arc = Arc(tour[i], tour[i+1])
        haskey(arc_idx, a) || continue
        kx::Int = arc_idx[a]
        set_start_value(w[kx], min(cum, Tlim))
        cum += arc_time(data, a)
        v::Int = tour[i+1]
        for (b::Int, vv::Int) in y_on
            vv == v && (cum += service_time(data, b))
        end
        cum = min(cum, Tlim)
    end
    return nothing
end

"""Parse optional `app[\"path_cbrp_fix_warm_start\"]` — Boolean or `\"true\"` / `\"false\"` strings."""
function pathCbrpFixWarmStartFlag(app::Dict{String,Any})::Bool
    v = get(app, "path_cbrp_fix_warm_start", false)
    v isa Bool && return v
    v isa AbstractString && return lowercase(strip(String(v))) in ("1", "true", "yes")
    try
        return Bool(v)
    catch
        return false
    end
end

"""
    pathCbrpFixWarmFromXY!(model, A, y_meta, x_on, y_on)::Nothing

Hard-fix Path-CBRP binaries `x` and `y` with explicit equalities `fix_warm_x[k]` / `fix_warm_y[k]` matching sparse arc set `x_on` and assignments `y_on` (solver-visible names for diagnostics).
Does not fix `w`; solver must find feasible MTZ potentials if they exist.
"""
function pathCbrpFixWarmFromXY!(
    model::Model,
    A::Arcs,
    y_meta::Vector{Tuple{Int,Int}},
    x_on::Set{Arc},
    y_on::Set{Tuple{Int,Int}},
)::Nothing
    na::Int = length(A)
    ny::Int = length(y_meta)
    arc_idx::Dict{Arc,Int} = Dict{Arc,Int}(a => k for (k, a) in enumerate(A))
    x = model[:x]
    y = model[:y]

    for a::Arc in x_on
        haskey(arc_idx, a) || throw(ArgumentError("Fix warm-start arc $a not in sparse arc list"))
    end
    y_meta_set::Set{Tuple{Int,Int}} = Set(y_meta)
    for pair::Tuple{Int,Int} in y_on
        pair in y_meta_set || throw(ArgumentError("Fix warm-start y assignment $pair not in y_meta"))
    end

    rhs_x::Vector{Float64} = [Float64(A[k] ∈ x_on) for k in 1:na]
    rhs_y::Vector{Float64} = [Float64(y_meta[k] ∈ y_on) for k in 1:ny]
    @constraint(model, fix_warm_x[k=1:na], x[k] == rhs_x[k])
    @constraint(model, fix_warm_y[k=1:ny], y[k] == rhs_y[k])
    return nothing
end

"""
    pathCbrpFixWarmFromSolution!(model, data, A, y_meta, sol)::Nothing

Fix `x`/`y` via explicit equalities `fix_warm_x[k]` / `fix_warm_y[k]` from a sparse `SBRPSolution` tour and serviced blocks (same semantics as `pathCbrpMipStart!`).
"""
function pathCbrpFixWarmFromSolution!(
    model::Model,
    data::SBRPData,
    A::Arcs,
    y_meta::Vector{Tuple{Int,Int}},
    sol::SBRPSolution,
)::Nothing
    na::Int = length(A)
    ny::Int = length(y_meta)
    arc_idx::Dict{Arc,Int} = Dict{Arc,Int}(a => k for (k, a) in enumerate(A))
    x = model[:x]
    y = model[:y]
    tour::Vi = sol.tour

    x_used::Set{Arc} = Set{Arc}()
    for i::Int in 1:(length(tour)-1)
        a::Arc = Arc(tour[i], tour[i+1])
        haskey(arc_idx, a) || throw(ArgumentError("Fix warm-start tour arc $a not in sparse arc list"))
        push!(x_used, a)
    end

    served_y::Set{Tuple{Int,Int}} = Set{Tuple{Int,Int}}()
    for b::Int in 1:num_clusters(data)
        blk::Vi = get_cluster(data, b)
        any(cand::Vi -> Set(cand) == Set(blk), sol.B) || continue
        v_pick::Union{Nothing,Int} = nothing
        for node::Int in tour
            if node in blk
                v_pick = node
                break
            end
        end
        v_pick === nothing &&
            throw(ArgumentError("Fix warm start: served block $b has no tour vertex in that block"))
        push!(served_y, (b, v_pick))
    end

    rhs_x::Vector{Float64} = [Float64(A[k] ∈ x_used) for k in 1:na]
    rhs_y::Vector{Float64} = [Float64(y_meta[k] ∈ served_y) for k in 1:ny]
    @constraint(model, fix_warm_x[k=1:na], x[k] == rhs_x[k])
    @constraint(model, fix_warm_y[k=1:ny], y[k] == rhs_y[k])
    return nothing
end

"""
    detectPathCbrpWarmSolFormat(path)::String

Return `path-cbrp-mtz` if the file contains `X:` lines, else `brkga-sol`.
"""
function detectPathCbrpWarmSolFormat(path::String)::String
    for line::String in readlines(path)
        startswith(strip(line), "X:") && return "path-cbrp-mtz"
    end
    return "brkga-sol"
end

"""
    resolvePathCbrpWarmSolFormat(app, path)::String

Resolve `brkga-sol`, `path-cbrp-mtz`, or `auto` (default).
"""
function resolvePathCbrpWarmSolFormat(app::Dict{String,Any}, path::String)::String
    fmt::String = String(strip(String(get(app, "path-cbrp-warm-sol-format", "auto"))))
    fmt == "auto" && return detectPathCbrpWarmSolFormat(path)
    fmt in ("brkga-sol", "path-cbrp-mtz") ||
        error("Unknown --path-cbrp-warm-sol-format: $(repr(fmt)) (use brkga-sol, path-cbrp-mtz, or auto)")
    return fmt
end

"""
    parsePathCbrpMtzSolution(path)::Tuple{Set{Arc}, Vector{Tuple{Int,Int}}}

Parse article `path-cbrp-mtz` files: `X: i j` arcs and `Y: node block_id` assignments.
"""
function parsePathCbrpMtzSolution(path::String)::Tuple{Set{Arc},Vector{Tuple{Int,Int}}}
    x_on::Set{Arc} = Set{Arc}()
    y_node_block::Vector{Tuple{Int,Int}} = Tuple{Int,Int}[]
    for line::String in readlines(path)
        s::String = strip(line)
        isempty(s) && continue
        if startswith(s, "X:")
            parts::Vector{String} = split(s, [' ', ':', '\t']; keepempty=false)
            length(parts) >= 3 || throw(ArgumentError("Invalid X line: $line"))
            push!(x_on, Arc(parse(Int, parts[2]), parse(Int, parts[3])))
        elseif startswith(s, "Y:")
            parts = split(s, [' ', ':', '\t']; keepempty=false)
            length(parts) >= 3 || throw(ArgumentError("Invalid Y line: $line"))
            push!(y_node_block, (parse(Int, parts[2]), parse(Int, parts[3])))
        end
    end
    return x_on, y_node_block
end

"""
    buildCarlosBlockIdToIndex(data, instance_path)::Dict{Int,Int}

Map Carlos instance file `block_id` to Julia cluster index `b` (1..num_clusters)
by matching block vertex sets to `data.B[b]`.
"""
function buildCarlosBlockIdToIndex(data::SBRPData, instance_path::String)::Dict{Int,Int}
    block_arcs::Dict{Int,Arcs} = Dict{Int,Arcs}()
    open(instance_path) do f::IOStream
        hdr::Vector{Int} = map(
            str::String -> parse(Int, str),
            Vector{String}(split(readline(f), [' ']; limit=3, keepempty=false)),
        )
        nNodes::Int = hdr[1]
        nArcs::Int = hdr[2]
        for _ in 1:nNodes
            readline(f)
        end
        for _ in 1:nArcs
            strs::Vector{String} = split(readline(f), [' ']; keepempty=false)
            block_id::Int = parse(Int, strs[5])
            block_id == -1 && continue
            a::Arc = Arc(parse(Int, strs[2]), parse(Int, strs[3]))
            if !haskey(block_arcs, block_id)
                block_arcs[block_id] = Arcs()
            end
            push!(block_arcs[block_id], a)
        end
    end

    id_to_b::Dict{Int,Int} = Dict{Int,Int}()
    for (block_id::Int, arcs::Arcs) in block_arcs
        arcs_mut::Arcs = Arcs(copy(arcs))
        block::Vi = _carlosBlockVerticesFromArcs!(arcs_mut)
        b::Union{Int,Nothing} = findfirst(
            bb::Int -> Set{Int}(get_cluster(data, bb)) == Set{Int}(block),
            1:num_clusters(data),
        )
        b === nothing &&
            throw(ArgumentError("Carlos block_id $block_id not found in loaded instance data.B"))
        id_to_b[block_id] = b
    end
    return id_to_b
end

function _attachPathCbrpWarmStartBrkgaSol!(app::Dict{String,Any}, data_sparse::SBRPData, path_str::String)::Nothing
    app_dense::Dict{String,Any} = copy(app)
    app_dense["no-cbrp-metric-closure"] = false
    data_compact::SBRPData, paths_u, dist_u = readSBRPDataCarlos(app_dense)
    paths_u === nothing && error("Path-CBRP warm start: metric closure did not yield paths")
    dist_u === nothing && error("Path-CBRP warm start: metric closure did not yield street distances")
    paths::Dict{Arc,Vi} = paths_u
    street_distances::ArcCostMap = dist_u
    data_compact.T = data_sparse.T
    sol_compact::SBRPSolution = readSBRPSolution(data_compact, path_str, street_distances)
    app["warm_start_solution"] = compactToSparseSbrpSolution(data_sparse, paths, sol_compact)
    return nothing
end

function _attachPathCbrpWarmStartMtz!(app::Dict{String,Any}, data_sparse::SBRPData, path_str::String)::Nothing
    x_on::Set{Arc}, y_nb::Vector{Tuple{Int,Int}} = parsePathCbrpMtzSolution(path_str)
    instance_path::String = String(app["instance"])
    id_to_b::Dict{Int,Int} = buildCarlosBlockIdToIndex(data_sparse, instance_path)

    y_on::Set{Tuple{Int,Int}} = Set{Tuple{Int,Int}}()
    seen_block::Dict{Int,Int} = Dict{Int,Int}()
    for (node::Int, block_id::Int) in y_nb
        haskey(id_to_b, block_id) ||
            throw(ArgumentError("Y references unknown Carlos block_id $block_id"))
        b::Int = id_to_b[block_id]
        blk::Vi = get_cluster(data_sparse, b)
        node in blk || throw(ArgumentError("Y: node $node not in block $block_id (Julia b=$b)"))
        if haskey(seen_block, b) && seen_block[b] != node
            @warn "Path-CBRP warm start: multiple Y for Julia block $b (nodes $(seen_block[b]) and $node); keeping first"
            continue
        end
        seen_block[b] = node
        push!(y_on, (b, node))
    end

    for a::Arc in x_on
        haskey(data_sparse.D.distance, a) ||
            throw(ArgumentError("X arc $a not in sparse instance digraph"))
    end

    app["path_cbrp_warm_xy"] = (x_on, y_on)
    return nothing
end

"""
    attachPathCbrpWarmStartFromFile!(app, data_sparse)::Nothing

Load Path-CBRP MIP warm start from `app[\"path-cbrp-warm-sol\"]` using
`app[\"path-cbrp-warm-sol-format\"]` (`brkga-sol`, `path-cbrp-mtz`, or `auto`).
"""
function attachPathCbrpWarmStartFromFile!(app::Dict{String,Any}, data_sparse::SBRPData)::Nothing
    path_str::String = String(strip(String(get(app, "path-cbrp-warm-sol", ""))))
    isempty(path_str) && return nothing
    isfile(path_str) || error("Path-CBRP warm start: file not found: $(repr(path_str))")

    fmt::String = resolvePathCbrpWarmSolFormat(app, path_str)
    try
        if fmt == "brkga-sol"
            _attachPathCbrpWarmStartBrkgaSol!(app, data_sparse, path_str)
        else
            _attachPathCbrpWarmStartMtz!(app, data_sparse, path_str)
        end
    catch e
        @warn "Path-CBRP warm start failed ($fmt); continuing without MIP start: $(sprint(showerror, e))"
        haskey(app, "warm_start_solution") && delete!(app, "warm_start_solution")
        haskey(app, "path_cbrp_warm_xy") && delete!(app, "path_cbrp_warm_xy")
    end
    return nothing
end

"""Deprecated alias."""
attachPathCbrpWarmStartFromBrkgaSolFile!(app::Dict{String,Any}, data_sparse::SBRPData) =
    attachPathCbrpWarmStartFromFile!(app, data_sparse)
