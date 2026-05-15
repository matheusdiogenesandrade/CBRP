include("symbols.jl")
include("data.jl")
include("sol.jl")
include("nn_heuristic.jl")
include("model.jl")
include("brkga.jl")

using ArgParse
using Logging

# log
const LOG_FILE = open("logs/log", "w+")
global_logger(ConsoleLogger(LOG_FILE, Logging.Debug))
disable_logging(Debug)

function parse_commandline(args_array::Vector{String}, appfolder::String)::Union{Nothing,Dict{String,Any}}
    s::ArgParseSettings = ArgParseSettings(usage="  On interactive mode, call main([\"arg1\", ..., \"argn\"])", exit_after_help=false)

    @add_arg_table s begin
        "instance"
        help = "Instance file path"
        "--cluster-size-profits"
        help = "true if you want to consider the blocks' profits as the size of the block, and false otherwise"
        action = :store_true
        "--unitary-profits"
        help = "true if you want to consider the blocks' profits as 1, and false otherwise"
        action = :store_true
        "--sol-stats"
        help = "true if you want to generate solution stats from an input solution, and false otherwise"
        action = :store_true
        "--ip"
        help = "true if you want to run the I.P. model, and false otherwise"
        action = :store_true
        "--brkga"
        help = "true if you want to run the BRKGA, and false otherwise"
        action = :store_true
        "--brkga-conf"
        help = "BRKGA config file directory"
        default = "conf/config.conf"
        "--vehicle-time-limit"
        help = "Vehicle time limit in minutes"
        default = "120"
        "--time-limit"
        help = "Seconds budget for Path-CBRP CPLEX solve (0 → 3600). Other IP modes unchanged."
        default = "0"
        "--instance-type"
        help = "Instance type (matheus|carlos)"
        default = "matheus"
        "--nosolve"
        help = "Not solve flag"
        action = :store_true
        "--out"
        help = "Path to write the solution found"
        "--batch"
        help = "Batch file path"
        "--intersection-cuts"
        help = "Intersection cuts for the complete model"
        action = :store_true
        "--subcycle-separation"
        help = "Strategy (first|best|all|none)"
        help = "Subcycle separation with max-flow at the B&B root node"
        default = "all"
        "--y-integer"
        help = "Fix the variable y, for the complete model, when running the separation algorithm"
        action = :store_true
        "--z-integer"
        help = "Fix the variable z, for the complete model, when running the separation algorithm"
        action = :store_true
        "--w-integer"
        help = "Fix the variable w, for the complete model, when running the separation algorithm"
        action = :store_true
        "--no-cbrp-metric-closure"
        help = "Carlos: skip Floyd–Warshall metric closure; keep sparse street digraph (pair with Path-CBRP IP)"
        action = :store_true
        "--drop-zero-profit-blocks"
        help = "Remove blocks with profit exactly 0 before metric closure / complete digraph (Carlos and Matheus readers)"
        action = :store_true
        "--path-cbrp-mip"
        help = "Carlos + IP: arc-indexed Path-CBRP MILP (requires --no-cbrp-metric-closure)"
        action = :store_true
        "--path-cbrp-warm-sol"
        help = "Carlos Path-CBRP: warm-start solution file (BRKGA .sol or path-cbrp-mtz article .txt)"
        arg_type = String
        default = ""
        "--path-cbrp-warm-sol-format"
        help = "Warm-start file format: brkga-sol | path-cbrp-mtz | auto (default auto: detect X: lines)"
        arg_type = String
        default = "auto"
    end

    return parse_args(args_array, s)
end

# log function
function log(app::Dict{String,Any}, info::Dict{String,String})

    columns::Vector{String} = ["instance", "|V|", "|A|", "|B|", "T", "model", "initialLP", "yLP", "yLPTime", "zLP", "zLPTime", "wLP", "wLPTime", "maxFlowLP", "maxFlowCuts", "maxFlowCutsTime", "lazyCuts", "cost", "solverTime", "relativeGAP", "nodeCount", "integerCount", "phase1Time", "meters", "tourMinutes", "blocksMeters", "blocksMinutes", "numVisitedBlocks", "intersectionCutsTime", "intersectionCuts1", "intersectionCuts2", "numVisitedNodes", "numOriginalVisitedNodes", "numRepeatedNodes", "numRepeatedArcs", "avgDetourIndex"]

    info["instance"] = last(split(app["instance"], "/"; keepempty=false))
    info["instance"] = first(split(info["instance"], "."; keepempty=false))

    selected_columns = filter(column -> column in keys(info), columns)

    # log
    @info join(selected_columns, ",")
    @info join(map(column -> info[column], selected_columns), ",")

    # console dump
    println(join(selected_columns, ","))
    println(join(map(column -> info[column], selected_columns), ","))
    flush(stdout)

end

#=
Get solution for original graph
    input:
        - tour::Vector{Integer} is the solution route
        - A::Vector{Pair{Int, Int}} is the list of arcs `tour` was built on
        - paths::Dict{Pair{Int, Int}, Vector{Int}} is the relation of paths, from the original digraph, between two nodes belonging to the blocks

    output:
        - ori_tour::Vector{Int} is the solution route using arcs from the original digraph
=#
function retrieveOriginalDigraphSolution(
    tour::Vi,
    paths::Dict{Arc,Vi},
    distances::Dict{Arc,Float64}
)::Vi

    # original route
    ori_tour::Vi = Vi()  # start
    total_distance::Float64 = 0.0
    for i in 1:length(tour)-1
        j = i + 1
        arc::Arc = Arc(tour[i], tour[j])
        path::Vi = paths[arc]
        #	append!(ori_tour, path[1:end-1])  # avoid duplication of nodes
        append!(ori_tour, tour[i])
        append!(ori_tour, path)
    end

    push!(ori_tour, tour[end])

    # return
    return ori_tour
end

function getPathsForIntersectionCuts(data::SBRPData, info::Dict{String,String})::Tuple{VVi,ArcsSet}

    # get
    elapsed_time::Float64 = @elapsed maximal_paths::VVi, invalid_arcs::ArcsSet = getMaximalPathsAndInvalidArcs(data)

    # paths statistics
    paths_lengths::Vi = map(maximal_path::Vi -> length(maximal_path), maximal_paths)

    num_paths::Int = length(maximal_paths)
    average::Float64 = reduce(+, paths_lengths) / num_paths
    max_length::Int = maximum(paths_lengths)
    min_length::Int = minimum(paths_lengths)
    std::Float64 = reduce(+, map(path_length::Int -> abs(path_length - average), paths_lengths)) / num_paths

    # store
    info["#PreprocArcs"] = string(length(invalid_arcs))
    info["#IntersectPaths"] = string(num_paths)
    info["AVGIntersectLen"] = string(average)
    info["MaxIntersectLen"] = string(max_length)
    info["MinIntersectLen"] = string(min_length)
    info["STDIntersectLen"] = string(std)
    info["IntersecTime"] = string(elapsed_time)

    # return
    return maximal_paths, invalid_arcs
end

#=
Run complete digraph IP formulation for the SBRP
    input:
        - app::Dict{String, Any} is the command line arguments relation
        - data::SBRPData is the SBRP instance built on a complete digraph
=#
function completeDigraphIPModel(
    app::Dict{String,Any},
    data::SBRPData,
    paths::Union{Dict{Arc,Vi},Nothing},
    distances::Union{ArcCostMap,Nothing},
)

    @info "###################SBRP####################"

    # intersection cuts
    info::Dict{String,String} = Dict{String,String}()
    maximal_paths::VVi = VVi()

    # create and solve model
    solution::SBRPSolution, info_::Dict{String,String} = runCOPCompleteDigraphIPModel(data, app)

    # merge
    merge!(info, info_)

    # check feasibility
    checkSBRPSolution(data, solution)

    # log
    info["model"] = "IP"
    info["|V|"] = string(length(data.D.V))
    info["|A|"] = string(length(data.D.A))
    info["|B|"] = string(length(data.B))
    info["T"] = string(data.T)
    info["meters"] = string(tourDistance(data, solution.tour))
    info["tourMinutes"] = string(tourTime(data, solution))
    info["blocksMeters"] = string(sum(map(block::Vi -> blockDistance(data, block), solution.B)))
    info["numVisitedBlocks"] = string(length(solution.B))

    log(app, info)

    # write solution
    solution_dir::Union{String,Nothing} = app["out"]

    if solution_dir != nothing
        if paths !== nothing && distances !== nothing
            solution.tour = retrieveOriginalDigraphSolution(solution.tour[2:end-1], paths, distances)
        else
            @warn "Skipping compact→street tour expansion (--out): paths/street distances not available (e.g. Carlos sparse reader)"
        end

        writeSolution(solution_dir * "_cop_ip_model", data, solution)
    end

    @info "########################################################"
end

#=
Path-CBRP IP (Carlos sparse digraph): arc-indexed MILP, no complete-graph MTZ.
=#
function pathCbrpIPModel(app::Dict{String,Any}, data::SBRPData)

    @info "###################Path-CBRP (Carlos sparse)####################"

    app["path_cbrp_fix_warm_start"] = true

    solution::SBRPSolution, info_::Dict{String,String} = runPathCbrpMipModel(data, app)
    info::Dict{String,String} = Dict{String,String}()
    merge!(info, info_)

    try
        checkSBRPSolution(data, solution)
    catch e
        @warn "Path-CBRP solution failed SBRP arc/block checks (subtours possible without SEC): $(sprint(showerror, e))"
    end

    info["model"] = "PathCBRPIP"
    info["|V|"] = string(length(data.D.V))
    info["|A|"] = string(length(data.D.A))
    info["|B|"] = string(length(data.B))
    info["T"] = string(data.T)
    info["meters"] = string(tourDistance(data, solution.tour))
    info["tourMinutes"] = string(tourTime(data, solution))
    info["blocksMeters"] = string(sum(map(block::Vi -> blockDistance(data, block), solution.B); init=0.0))
    info["numVisitedBlocks"] = string(length(solution.B))

    log(app, info)

    solution_dir::Union{String,Nothing} = app["out"]
    if solution_dir != nothing
        writeSolution(solution_dir * "_cop_ip_model", data, solution)
    end

    @info "########################################################"
end

#=
Run BRKGA algorithm
    input:
        - app::Dict{String, Any} is the command line arguments relation
        - data::SBRPData is the SBRP instance built on a complete digraph
=#
function BRKGAModel(app::Dict{String,Any}, data::SBRPData)

    @info "###################BRKGA####################"

    # solve model
    solution::SBRPSolution, info::Dict{String,String} = runCOPBRKGAModel(data, app)

    # check feasibility
    checkSBRPSolution(data, solution)

    # log
    info["model"] = "BRKGA"
    info["|V|"] = string(length(data.D.V))
    info["|A|"] = string(length(data.D.A))
    info["|B|"] = string(length(data.B))
    info["T"] = string(data.T)
    info["meters"] = string(tourDistance(data, solution.tour))
    info["tourMinutes"] = string(tourTime(data, solution))
    info["blocksMeters"] = string(sum(map(block::Vi -> blockDistance(data, block), solution.B)))
    info["numVisitedBlocks"] = string(length(solution.B))

    log(app, info)

    # write solution
    solution_dir::Union{String,Nothing} = app["out"]
    if solution_dir != nothing
        writeSolution(solution_dir * "_brkga", data, solution)
    end

    @info "########################################################"
end

#=
Calculate the number of times a tour crosses between different communities.
=#
function calculate_inter_community_crossings(tour::Vi, communities::Dict{Int,Int})::Int
    crossings = 0
    if length(tour) > 1
        for i in 1:(length(tour)-1)
            u, v = tour[i], tour[i+1]
            if communities[u] != communities[v]
                crossings += 1
            end
        end
    end
    return crossings
end

#=
Fill solution statistics
input:
    - app::Dict{String, Any} is the command line arguments relation
    - data::SBRPData is the SBRP instance built on a complete digraph
    - solution::SBRPSolution is the solution to be filled
    - paths::Dict{Arc, Vi} is the relation of paths, from the original digraph, between two nodes belonging to the blocks
    - distances::Dict{Arc, Float64} is the distance between two nodes in the original digraph
=#
function fillSolutionStats(app::Dict{String,Any}, data::SBRPData, paths::Union{Dict{Arc,Vi},Nothing}, distances::Union{Dict{Arc,Float64},Nothing})

    @info "###################Solution-Statistics####################"

    # check paths
    if paths == nothing
        error("You must provide the paths when using the --sol-stats option")
    end

    # check distances
    if distances == nothing
        error("You must provide the distances when using the --sol-stats option")
    end

    # read solution
    solution_dir::Union{String,Nothing} = app["out"]

    # check
    if solution_dir == nothing
        error("You must provide a solution file path using the --out option")
    end

    # info
    info::Dict{String,String} = Dict{String,String}()

    # read
    solution::SBRPSolution = readSBRPSolution(data, solution_dir, distances)

    # log
    info["model"] = "SOL-STATS"
    info["|V|"] = string(length(data.D.V))
    info["|A|"] = string(length(data.D.A))
    info["|B|"] = string(length(data.B))
    info["T"] = string(data.T)

    info["cost"] = string(sum(block::Vi -> data.profits[block], solution.B; init=0.0))


    info["totalBlocksNodes"] = string(sum(block::Vi -> length(block), solution.B; init=0.0))

    info["meters"] = string(tourDistance(distances, solution.tour))
    info["blocksMeters"] = string(sum(map(block::Vi -> blockDistance(data, distances, block), solution.B); init=0.0))

    info["blocksMinutes"] = string(sum(block::Vi -> blockTime(data, distances, block), solution.B; init=0.0))
    info["tourMinutes"] = string(tourTime(data, distances, solution) - sum(block::Vi -> blockTime(data, distances, block), solution.B; init=0.0))

    info["numVisitedBlocks"] = string(length(solution.B))

    #
    info["numVisitedNodes"] = string(count_nodes_in_serviced_blocks(solution.tour, solution.B))
    info["numOriginalVisitedNodes"] = string(length(solution.tour))

    #
    info["numRepeatedNodes"] = string(getNumRepeatedNodes(solution.tour))
    info["numRepeatedArcs"] = string(getNumRepeatedArcs(solution.tour))
    #    info["avgDetourIndex"] = string(calculate_average_detour_index(solution.tour, solution, data.D.distance))

    #
    log(app, info)

    @info "########################################################"
end

function run(app::Dict{String,Any})
    @info "Application parameters:"

    for (arg, val) in app
        @info "  $arg  =>  $(repr(val))"
    end

    if app["path-cbrp-mip"]
        app["instance-type"] == "carlos" || error("--path-cbrp-mip requires --instance-type carlos")
        app["ip"] || error("--path-cbrp-mip requires --ip")
        app["no-cbrp-metric-closure"] || error("--path-cbrp-mip requires --no-cbrp-metric-closure")
        app["brkga"] && error("--path-cbrp-mip cannot be combined with --brkga")
    end

    warm_sol_path::String = String(strip(String(get(app, "path-cbrp-warm-sol", ""))))
    if !isempty(warm_sol_path)
        app["path-cbrp-mip"] || error("--path-cbrp-warm-sol requires --path-cbrp-mip")
        app["ip"] || error("--path-cbrp-warm-sol requires --ip")
        app["instance-type"] == "carlos" || error("--path-cbrp-warm-sol requires --instance-type carlos")
        app["no-cbrp-metric-closure"] || error("--path-cbrp-warm-sol requires --no-cbrp-metric-closure")
        app["brkga"] && error("--path-cbrp-warm-sol cannot be combined with --brkga")
        fmt_chk::String = String(strip(String(get(app, "path-cbrp-warm-sol-format", "auto"))))
        fmt_chk in ("brkga-sol", "path-cbrp-mtz", "auto") ||
            error("Unknown --path-cbrp-warm-sol-format: $(repr(fmt_chk))")
    end

    # read instance
    data::Union{SBRPData,Nothing} = nothing
    paths::Union{Dict{Arc,Vi},Nothing} = nothing
    distances::Union{Dict{Arc,Float64},Nothing} = nothing

    if app["instance-type"] == "matheus"
        data, paths, distances = readSBRPData(app)
    elseif app["instance-type"] == "carlos"
        data, paths, distances = readSBRPDataCarlos(app)
    else
        error(@sprintf("Invalid instance type: %s", app["instance-type"]))
    end

    # instance data
    @info "|B| = $(length(data.B))"
    @info "|V| = $(length(data.D.V))"
    @info "|A| = $(length(data.D.A))"

    # set vehicle time limit
    data.T = parse(Int, app["vehicle-time-limit"])

    # not solve
    if app["nosolve"]

        info::Dict{String,String} = Dict{String,String}()

        # log
        info["|V|"] = string(length(data.D.V))
        info["|A|"] = string(length(data.D.A))
        info["|B|"] = string(length(data.B))
        info["T"] = string(data.T)

        log(app, info)
    elseif app["ip"]
        if app["path-cbrp-mip"]
            attachPathCbrpWarmStartFromFile!(app, data)
            pathCbrpIPModel(app, data)
            haskey(app, "warm_start_solution") && delete!(app, "warm_start_solution")
            haskey(app, "path_cbrp_warm_xy") && delete!(app, "path_cbrp_warm_xy")
        else
            completeDigraphIPModel(app, data, paths, distances)
        end
    elseif app["brkga"]
        BRKGAModel(app, data)
    elseif app["sol-stats"]
        fillSolutionStats(app, data, paths, distances)
    else
        error("You must choose at least one of: --ip, --brkga, or --sol-stats")
    end
end

function main(args)
    appfolder::String = dirname(@__FILE__)

    app::Union{Nothing,Dict{String,Any}} = parse_commandline(args, appfolder)

    isnothing(app) && return

    if app["batch"] != nothing
        for line in readlines(app["batch"])
            if !isempty(strip(line)) && strip(line)[1] != '#'
                #		    try
                run(parse_commandline(map(s::Any -> String(s), split(line)), appfolder))
                #		    catch e
                #			    @error "Error processing line: $line"
                #			    @error "$e"
                #			    @error "$(catch_backtrace())"
                #		    end
                GC.gc()  # Force the GC to run

            end
        end
    else
        run(app)
    end
end

if isempty(ARGS)
    main(["--help"])
else
    main(ARGS)
end
