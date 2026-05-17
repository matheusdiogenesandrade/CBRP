#!/usr/bin/env julia
#
# Complete digraph IP time-budget sweep (Carlos + metric closure): repeatedly solve
# runCOPCompleteDigraphIPModel while increasing data.T until num_serviced_blocks >=
# num_positive_profit_blocks, a solver exception occurs, or caps (max-T, max-iterations) apply.
#
# Requires Carlos with metric closure (omit sparse flag): uses the same reader as
# `src/run.jl` with default closure for complete digraph IP.
#
# Cut / warm-start pooling (optional):
#   --reuse-cuts accumulates subtour (max-flow) cuts across budgets, caches intersection
#   cuts after the first run (when --intersection-cuts), and warm-starts the MIP from the
#   best prior feasible solution (single JuMP primal hint). Max-flow separation always runs
#   when subcycle separation is enabled (pooled cuts are seeded first, then separation continues).
#
# Usage (from external/CBRP_original):
#   julia --threads=1 --project=. experiments/complete_digraph_time_sweep.jl INSTANCE.txt --out-csv OUT.csv
#
# Optional --drop-zero-profit-blocks: remove blocks with profit exactly 0 before metric
# closure (same as src/run.jl). CSV num_blocks / num_positive_profit_blocks refer to the
# reduced instance passed to the IP.

using ArgParse

const ROOT = normpath(joinpath(@__DIR__, ".."))

include(joinpath(ROOT, "src", "data.jl"))
include(joinpath(ROOT, "src", "model.jl"))

function csv_escape(s::AbstractString)::String
    if occursin(r"[\",\n\r]", s)
        return "\"" * replace(s, "\"" => "\"\"") * "\""
    end
    return String(s)
end

function parse_cli()::Union{Nothing,Dict{String,Any}}
    s = ArgParseSettings(
        prog = "complete_digraph_time_sweep.jl",
        description = "Sweep complete digraph IP (Carlos, metric closure) over increasing time budget T. Stops when num_serviced_blocks >= count of blocks with profit > 0 on the instance actually solved (after optional --drop-zero-profit-blocks).",
        exit_after_help = false,
    )
    @add_arg_table s begin
        "instance"
        help = "Carlos instance .txt path"
        arg_type = String
        "--min-T"
        help = "First budget T in minutes"
        arg_type = Float64
        default = 10.0
        "--time-step"
        help = "Increase T by this many minutes each iteration"
        arg_type = Float64
        default = 5.0
        "--max-T"
        help = "Inclusive upper cap on budget T (minutes); stop if next T would exceed this"
        arg_type = Float64
        default = 400.0
        "--max-iterations"
        help = "Maximum number of complete IP solves"
        arg_type = Int
        default = 100
        "--time-limit"
        help = "CPLEX per-solve time limit in seconds (0 → 3600, same as src/run.jl)"
        arg_type = String
        default = "0"
        "--subcycle-separation"
        help = "Subcycle separation: first|best|all|none (complete IP prep phase)"
        arg_type = String
        default = "all"
        "--intersection-cuts"
        help = "Enable intersection cuts in the complete model"
        action = :store_true
        "--y-integer"
        help = "Keep y integer during subtour separation (see src/run.jl)"
        action = :store_true
        "--z-integer"
        help = "Keep z integer during subtour separation"
        action = :store_true
        "--w-integer"
        help = "Keep w integer during subtour separation"
        action = :store_true
        "--reuse-cuts"
        help = "Accumulate subtour cuts, cache intersection cuts after first solve, warm-start from best prior solution"
        action = :store_true
        "--warm-start-pool-size"
        help = "Max number of prior feasible solutions kept for picking the best warm start (requires --reuse-cuts)"
        arg_type = Int
        default = 5
        "--out-csv"
        help = "Output CSV file path"
        arg_type = String
        required = true
        "--drop-zero-profit-blocks"
        help = "Remove blocks with profit exactly 0 before complete digraph (see src/sbrp.jl dropZeroProfitBlocks!)"
        action = :store_true
    end
    return parse_args(ARGS, s)
end

function get_info_str(info::Dict{String,String}, key::String)::String
    return get(info, key, "")
end

function write_result_row(
    io::IO,
    iteration::Int,
    budget_T::Float64,
    num_blocks::Int,
    num_positive_profit_blocks::Int,
    num_serviced::Int,
    terminal::String,
    info::Dict{String,String},
    error_message::String,
)
    cols = String[
        string(iteration),
        string(budget_T),
        string(num_blocks),
        string(num_positive_profit_blocks),
        string(num_serviced),
        csv_escape(terminal),
        csv_escape(get_info_str(info, "cost")),
        csv_escape(get_info_str(info, "bestBound")),
        csv_escape(get_info_str(info, "solverTime")),
        csv_escape(get_info_str(info, "relativeGAP")),
        csv_escape(get_info_str(info, "nodeCount")),
        csv_escape(get_info_str(info, "yLPTime")),
        csv_escape(get_info_str(info, "maxFlowCutsTime")),
        csv_escape(get_info_str(info, "maxFlowCuts")),
        csv_escape(get_info_str(info, "intersectionCutsTime")),
        csv_escape(get_info_str(info, "phase1Time")),
        csv_escape(get_info_str(info, "noFeasibleSolution")),
        csv_escape(error_message),
    ]
    println(io, join(cols, ","))
    flush(io)
end

"""Shallow copy of tour and blocks for storing in the solution pool."""
function copy_sbrp_solution(sol::SBRPSolution)::SBRPSolution
    return SBRPSolution(copy(sol.tour), deepcopy(sol.B))
end

function main()::Cint
    app = parse_cli()
    if app === nothing
        return Cint(0)
    end
    if !haskey(app, "instance") || isempty(strip(app["instance"]))
        println(stderr, "error: instance path is required")
        return Cint(2)
    end

    min_T::Float64 = app["min-T"]
    time_step::Float64 = app["time-step"]
    max_T::Float64 = app["max-T"]
    max_iterations::Int = app["max-iterations"]
    time_limit_str::String = app["time-limit"]
    out_csv::String = app["out-csv"]
    reuse_cuts::Bool = get(app, "reuse-cuts", false)
    warm_pool_cap::Int = app["warm-start-pool-size"]

    if warm_pool_cap < 1
        println(stderr, "error: --warm-start-pool-size must be >= 1")
        return Cint(2)
    end

    inst_raw::String = app["instance"]
    instance_path::String = isabspath(inst_raw) ? inst_raw : abspath(inst_raw)

    if min_T > max_T
        println(stderr, "error: min-T ($min_T) must be <= max-T ($max_T)")
        return Cint(2)
    end
    if time_step <= 0.0
        println(stderr, "error: time-step must be positive")
        return Cint(2)
    end
    if max_iterations < 1
        println(stderr, "error: max-iterations must be at least 1")
        return Cint(2)
    end

    subcycle::String = app["subcycle-separation"]
    if !(subcycle in ("first", "best", "all", "none"))
        println(stderr, "error: --subcycle-separation must be one of first|best|all|none")
        return Cint(2)
    end

    # Metric closure ON (default): do not set no-cbrp-metric-closure
    app_read = Dict{String,Any}(
        "instance" => instance_path,
        "vehicle-time-limit" => string(min_T),
    )
    if get(app, "drop-zero-profit-blocks", false)
        app_read["drop-zero-profit-blocks"] = true
    end
    data::SBRPData = readSBRPDataCarlos(app_read)[1]

    num_positive_profit_blocks::Int = count(b -> data.profits[b] > 0, data.B)
    if num_positive_profit_blocks == 0
        println(stderr, "error: no blocks with profit > 0; complete IP has nothing to optimize")
        return Cint(2)
    end

    subtour_cut_pool::Set{Tuple{Arcs, Arcs}} = Set{Tuple{Arcs, Arcs}}()
    intersection_cuts_cache_ref::Ref{Union{Nothing,Tuple{Vector{Arcs},Vector{Arcs}}}} =
        Ref{Union{Nothing,Tuple{Vector{Arcs},Vector{Arcs}}}}(nothing)
    solution_pool::Vector{Tuple{SBRPSolution, Float64}} = Tuple{SBRPSolution, Float64}[]
    warm_start_next::Union{Nothing, SBRPSolution} = nothing

    io = open(out_csv, "w")
    try
        header = [
            "iteration",
            "budget_T",
            "num_blocks",
            "num_positive_profit_blocks",
            "num_serviced_blocks",
            "terminal",
            "cost",
            "bestBound",
            "solverTime",
            "relativeGAP",
            "nodeCount",
            "yLPTime",
            "maxFlowCutsTime",
            "maxFlowCuts",
            "intersectionCutsTime",
            "phase1Time",
            "noFeasibleSolution",
            "error_message",
        ]
        println(io, join(header, ","))
        flush(io)

        curr_T::Float64 = min_T
        iter::Int = 0
        num_blocks::Int = length(data.B)

        while true
            iter += 1
            if curr_T > max_T + 1e-9
                println(stderr, "error: budget_T ($curr_T) exceeds max-T ($max_T)")
                return Cint(2)
            end

            data.T = curr_T

            solve_app::Dict{String, Any} = Dict{String, Any}(
                "time-limit" => time_limit_str,
                "subcycle-separation" => subcycle,
                "intersection-cuts" => get(app, "intersection-cuts", false),
                "y-integer" => get(app, "y-integer", false),
                "z-integer" => get(app, "z-integer", false),
                "w-integer" => get(app, "w-integer", false),
            )

            if reuse_cuts
                solve_app["reuse-cuts"] = true
                solve_app["subtour_cut_pool"] = subtour_cut_pool
                solve_app["intersection_cuts_cache_ref"] = intersection_cuts_cache_ref
                ic_cached = intersection_cuts_cache_ref[]
                if get(app, "intersection-cuts", false) && ic_cached !== nothing
                    solve_app["precomputed_intersection_cuts1"] = ic_cached[1]
                    solve_app["precomputed_intersection_cuts2"] = ic_cached[2]
                end
                if warm_start_next !== nothing
                    solve_app["warm_start_solution"] = warm_start_next
                end
            end

            si = try
                runCOPCompleteDigraphIPModel(data, solve_app)
            catch e
                msg = sprint(showerror, e)
                println(stderr, "Complete digraph IP sweep: solver error at iteration=$iter budget_T=$curr_T: $msg")
                empty_info = Dict{String, String}()
                write_result_row(
                    io,
                    iter,
                    curr_T,
                    num_blocks,
                    num_positive_profit_blocks,
                    0,
                    "error",
                    empty_info,
                    msg,
                )
                return Cint(1)
            end
            solution, info = si

            if reuse_cuts
                cost_v::Float64 = try
                    parse(Float64, info["cost"])
                catch
                    -Inf
                end
                push!(solution_pool, (copy_sbrp_solution(solution), cost_v))
                sort!(solution_pool, by = x -> x[2], rev = true)
                while length(solution_pool) > warm_pool_cap
                    pop!(solution_pool)
                end
                warm_start_next = length(solution_pool) > 0 ? copy_sbrp_solution(first(solution_pool)[1]) : nothing
            end

            n_serviced::Int = length(solution.B)

            term::String = ""
            if n_serviced >= num_positive_profit_blocks
                term = "all_positive_profit_blocks"
            elseif iter >= max_iterations
                term = "max_iterations"
            elseif curr_T + time_step > max_T
                term = "max_T"
            end

            write_result_row(io, iter, curr_T, num_blocks, num_positive_profit_blocks, n_serviced, term, info, "")

            if term != ""
                return Cint(0)
            end

            curr_T += time_step
        end
    finally
        close(io)
    end
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    exit(Int(main()))
end
