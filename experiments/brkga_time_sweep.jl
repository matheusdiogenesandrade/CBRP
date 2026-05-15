#!/usr/bin/env julia
#
# BRKGA time-budget sweep (Carlos + metric closure): repeatedly call runCOPBRKGAModel while
# increasing data.T until num_serviced_blocks >= num_positive_profit_blocks, an exception
# occurs, or caps (max-T, max-iterations) apply.
#
# Same CSV column layout as complete_digraph_time_sweep.jl for easy comparison; IP-only
# columns (yLPTime, maxFlowCuts*, intersectionCutsTime, relativeGAP, nodeCount, …) stay
# empty for BRKGA unless runCOPBRKGAModel fills them later.
#
# Usage (from external/CBRP_original):
#   julia --threads=1 --project=. experiments/brkga_time_sweep.jl INSTANCE.txt \\
#       --brkga-conf path/to/wrapper.conf --out-csv OUT.csv
#
# Optional: --time-limit SECONDS (string/float as for run.jl). If > 0, overrides BRKGA
# maximum_time from the wrapper conf file for each budget solve (same behavior as
# external/CBRP/src/brkga.jl). If 0, conf file maximum_time applies.
#
# Optional --drop-zero-profit-blocks: same as src/run.jl / complete_digraph_time_sweep.jl.

using ArgParse

const ROOT = normpath(joinpath(@__DIR__, ".."))

include(joinpath(ROOT, "src", "data.jl"))
include(joinpath(ROOT, "src", "nn_heuristic.jl"))
include(joinpath(ROOT, "src", "brkga.jl"))

function csv_escape(s::AbstractString)::String
    if occursin(r"[\",\n\r]", s)
        return "\"" * replace(s, "\"" => "\"\"") * "\""
    end
    return String(s)
end

function parse_cli()::Union{Nothing,Dict{String,Any}}
    s = ArgParseSettings(
        prog = "brkga_time_sweep.jl",
        description = "Sweep BRKGA (Carlos, metric closure) over increasing time budget T. Stops when num_serviced_blocks >= count of blocks with profit > 0 on the instance actually solved (after optional --drop-zero-profit-blocks).",
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
        help = "Maximum number of BRKGA solves"
        arg_type = Int
        default = 100
        "--time-limit"
        help = "BRKGA wall time per solve in seconds (0 → use maximum_time from --brkga-conf wrapper)"
        arg_type = String
        default = "0"
        "--brkga-conf"
        help = "BRKGA wrapper .conf path (same as src/run.jl --brkga-conf)"
        arg_type = String
        required = true
        "--out-csv"
        help = "Output CSV file path"
        arg_type = String
        required = true
        "--drop-zero-profit-blocks"
        help = "Remove blocks with profit exactly 0 before the sweep (see src/sbrp.jl dropZeroProfitBlocks!)"
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
    brkga_conf_raw::String = app["brkga-conf"]

    inst_raw::String = app["instance"]
    instance_path::String = isabspath(inst_raw) ? inst_raw : abspath(inst_raw)
    brkga_conf_path::String = isabspath(brkga_conf_raw) ? brkga_conf_raw : abspath(brkga_conf_raw)

    if !isfile(brkga_conf_path)
        println(stderr, "error: --brkga-conf file not found: $brkga_conf_path")
        return Cint(2)
    end

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
        println(stderr, "error: no blocks with profit > 0; BRKGA has nothing to optimize")
        return Cint(2)
    end

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

            solve_app::Dict{String,Any} = Dict{String,Any}(
                "brkga-conf" => brkga_conf_path,
                "time-limit" => time_limit_str,
            )

            si = try
                runCOPBRKGAModel(data, solve_app)
            catch e
                msg = sprint(showerror, e)
                println(stderr, "BRKGA time sweep: error at iteration=$iter budget_T=$curr_T: $msg")
                empty_info = Dict{String,String}()
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
