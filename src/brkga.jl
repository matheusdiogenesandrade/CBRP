using BrkgaMpIpr
using ConfParser
using Dates
using Printf

"""
    @enum StopRule

Controls stop criteria. Stops either when:
- a given number of `GENERATIONS` is given;
- or a `TARGET` value is found;
- or no `IMPROVEMENT` is found in a given number of iterations.
"""
@enum StopRule begin
    GENERATIONS = 0
    TARGET = 1
    IMPROVEMENT = 2
end

"""
    function parse(::Type{StopRule}, value::String)::StopRule

Parse `value` into a `StopRule`.
"""
function parseRule(::Type{StopRule}, value::String)::StopRule

    @debug "Parsing GA rules"

    local_value::Char = uppercase(strip(value)[1])

    dict::Dict{Char, StopRule} = Dict{Char, StopRule}('G' => GENERATIONS, 'T' => TARGET, 'I' => IMPROVEMENT)

    local_value in keys(dict) && return dict[local_value]

    throw(ArgumentError("cannot parse $value as StopRule"))
end

function dijkstra(data::SBRPData, idxs_blocks::Vi)::Tuple{Float64, Dict{Tuple{Int, Int}, Float64}}

    @debug "Dijkstra algorithm"


    # attrs

    @debug "Getting params"

    m::Int = length(idxs_blocks)
    B::VVi         = data.B
    A::Arcs        = data.D.A
    T::Float64     = data.T
    first_idx::Int = first(idxs_blocks) 

    # dijkstra distance <<idx_block, node>, min path time> stores the min path time from the depot to the node `node' in the block at `idx_block'
    consumed_times::Dict{Tuple{Int, Int}, Float64} = Dict{Tuple{Int, Int}, Float64}(reduce(vcat,
                                                                                           map(idx_block::Int -> 
                                                                                               map(i::Int -> (idx_block, i)::Tuple{Int, Int} => T + 1, B[idx_block]),
                                                                                               idxs_blocks
                                                                                              )
                                                                                          ))

    # initial times

    @debug "Initial times"

    first_block_time::Float64 = blockTime(data, B[first_idx])
    for i::Int in B[first_idx]
        consumed_times[(first_idx, i)] = first_block_time
    end

    # BFS

    @debug "BFS"

    for i::Int in 1:m - 1

        # get indexes

        @debug "Get indexes"

        idx_curr_block::Int = idxs_blocks[i]
        idx_next_block::Int = idxs_blocks[i + 1] 

        # get curr and next blocks

        @debug "Get curr and next block"

        curr_block::Vi = B[idx_curr_block]
        next_block::Vi = B[idx_next_block]

        # calculate next block service time

        @debug "Get next block service time"

        next_block_time::Float64 = blockTime(data, next_block)

        # intersecting nodes

        @debug "Get intersecting nodes"

        intersection::Si = Si(intersect(curr_block, next_block))

        for j::Int in intersection
            consumed_times[(idx_next_block, j)] = next_block_time + consumed_times[(idx_curr_block, j)]
        end

        # non intersecting nodes

        @debug "Get non-intersecting nodes"

        for i::Int in curr_block
            for j::Int in next_block
                if Arc(i, j) in keys(data.D.distance)

                    consumed_times[(idx_next_block, j)] = min(consumed_times[(idx_next_block, j)], consumed_times[(idx_curr_block, i)] + time(data, Arc(i, j)) + next_block_time)

                end
            end
        end
    end

    @debug "Get last block index"

    pos_last = findlast(idx::Int -> any(i::Int -> consumed_times[(idx, i)] <= T, B[idx]), idxs_blocks)
    if pos_last === nothing
        return T + 1.0, consumed_times
    end
    last_block_idx::Int = idxs_blocks[pos_last]

    @debug "Return"

    candidates::Vi = filter(i::Int -> consumed_times[(last_block_idx, i)] <= T, B[last_block_idx])
    if isempty(candidates)
        return T + 1.0, consumed_times
    end

    return minimum(map(i::Int -> consumed_times[(last_block_idx, i)], candidates)), consumed_times
end

function getDijkstraRoute(data::SBRPData, consumed_times::Dict{Tuple{Int, Int}, Float64}, idxs_blocks::Vi)::Vi

    @debug "Get dijkstra route"

    # attrs

    @debug "Get params"

    m_pos = findlast(idx::Int -> any(i::Int -> consumed_times[(idx, i)] <= data.T, data.B[idx]), idxs_blocks)
    if m_pos === nothing || m_pos == 0
        return Vi([data.depot])
    end
    m::Int = m_pos

    tour::Vi       = Vi()
    B::VVi         = data.B
    A::Arcs        = data.D.A
    T::Float64     = data.T
    first_idx::Int = first(idxs_blocks)
    last_idx::Int  = idxs_blocks[m]

    # initialize tour

    @debug "Initialize tour"

    push!(tour, B[last_idx][findmin(i -> consumed_times[(last_idx, i)], B[last_idx])[2]])

    # BFS

    @debug "BFS"

    for i::Int in reverse(2:m)

        # get indexes

        @debug "Get curr and next block indexes"

        idx_curr_block::Int = idxs_blocks[i]
        idx_prev_block::Int = idxs_blocks[i - 1]

        # get curr and prev blocks

        @debug "Get curr and next block"

        curr_block::Vi = B[idx_curr_block]
        prev_block::Vi = B[idx_prev_block]

        # get intersection

        @debug "Get intersection"

        intersection::Vi = ∩(curr_block, prev_block)

        # get block time

        @debug "Get curr block time"

        curr_block_time::Float64 = blockTime(data, B[idx_curr_block])

        # if no intersecting nodes push backwards candidate

        @debug "No intersecting nodes edge-case"

        j::Int = last(tour)
        if !in(j, intersection)
            candidates::Vi = filter(
                                    i::Int -> consumed_times[(idx_prev_block, i)] + time(data, Arc(i, j)) + curr_block_time <= consumed_times[(idx_curr_block, j)], 
                                    prev_block
                                   )
            if !isempty(candidates)
                push!(tour, first(candidates))
            end
        end

    end

    return reverse(tour)

end

function getIdxsBlocks(chromosome::Array{Float64}, data::SBRPData)::Vi
    # attrs
    m::Int = length(data.B)

    # get all the alleles with weight > 0.5 and sort them in reversed order 
    permutation::Array{Tuple{Float64, Int64}} = Array{Tuple{Float64, Int64}}(undef, m)

    for (idx::Int, key::Float64) in enumerate(chromosome)
        permutation[idx] = Tuple{Float64, Int}((key, idx))
    end

    sort!(permutation, rev = true)

    # return indexes of the selected blocks
    return map((key, idx)::Tuple{Float64, Int} -> idx, permutation)
end

function decode!(chromosome::Array{Float64}, data::SBRPData, rewrite::Bool)::Float64

    # attrs
    B::VVi = data.B

    # get min tour time
    idxs_blocks::Vi = getIdxsBlocks(chromosome, data)
    tour_time::Float64, consumed_times::Dict{Tuple{Int, Int}, Float64} = dijkstra(data, idxs_blocks)

    last_feasible = findlast(idx::Int -> any(i::Int -> consumed_times[(idx, i)] <= data.T, data.B[idx]), idxs_blocks)
    if last_feasible === nothing || tour_time > data.T
        return -Inf
    end
    idxs_blocks = idxs_blocks[begin:last_feasible]

    return sum(idx::Int -> data.profits[B[idx]], idxs_blocks)
end

function runCOPBRKGAModel(data::SBRPData, app::Dict{String, Any})::Tuple{SBRPSolution, Dict{String, String}}

    @debug "Run BRKGA"

    # attrs
    @debug "Getting instance data"

    B::VVi = data.B
    m::Int = length(B)
    verbose::Bool, info::Dict{String, String} = false, Dict{String, String}()

    STARTING_TIME = datetime2unix(now())

    ########################################
    # Load configuration file and show basic info.
    ########################################

    @debug "Load configuration file and show basic info"

    conf::ConfParse              = ConfParse(app["brkga-conf"])
    parse_conf!(conf)

    println(retrieve(conf, "stop_rule"))
    seed::Int                          = parse(Int64, retrieve(conf, "seed"))
    stop_rule::StopRule                = parseRule(StopRule, retrieve(conf, "stop_rule"))
    stop_argument::Union{Float64, Int} = parse(stop_rule == TARGET ? Float64 : Int64, retrieve(conf, "stop_argument"))
    maximum_time::Float64              = parse(Float64, retrieve(conf, "maximum_time"))

    time_limit_override::Float64 = parse(Float64, get(app, "time-limit", "0"))
    if time_limit_override > 0.0
        maximum_time = time_limit_override
    end

    verbose                            = parse(Bool, retrieve(conf, "verbose"))
    perform_evolution::Bool            = parse(Bool, retrieve(conf, "evolution"))

    if maximum_time <= 0.0
        error("Maximum time must be larger than 0.0. Given $maximum_time.")
    end

    ########################################
    # Load config file and show basic info.
    ########################################
    
    @debug "Load config file and show basic info"

    brkga_params, control_params = load_configuration(retrieve(conf, "brkg_params"))

    if verbose
        println("""
                ------------------------------------------------------
                > Experiment started at $(Dates.now())
                > Configuration: $(app["brkga-conf"])
                > Algorithm Parameters:
                """)
        if !perform_evolution
            println(">    - Simple multi-start: on (no evolutionary operators)")
        end
        # log BRKGA parameters
        for field in fieldnames(BrkgaParams)
            println(">  - $field $(getfield(brkga_params, field))")
        end
        # log control parameters
        for field in fieldnames(ExternalControlParams)
            println(">  - $field $(getfield(control_params, field))")
        end
        println("""
                > Seed: $seed
                > Stop rule: $stop_rule
                > Stop argument: $stop_argument
                > Maximum time (s): $maximum_time
                > Number of parallel threads for decoding: $(Threads.nthreads())
                ------------------------------------------------------""")
    end

    ########################################
    # Adjust BRKGA parameters
    ########################################

    verbose && println("\n[$(Dates.Time(Dates.now()))] Generating initial tour...")

    # Generate a greedy solution to be used as warm start for BRKGA.
    profit, visited_blocks, tour = solveNearestNeighborhood(data)
    verbose && println("Initial profit: $profit")

    ########################################
    # Build the BRKGA data structures and initialize
    ########################################

    verbose && println("\n[$(Dates.Time(Dates.now()))] Building BRKGA data...")

    # Usually, it is a good idea to set the population size
    # proportional to the instance size.
    brkga_params.population_size = min(brkga_params.population_size, 10 * m)
    verbose && println("New population size: $(brkga_params.population_size)")

    # Chromosome size is the number of nodes.
    # Each chromosome represents a permutation of nodes.
    brkga_data = build_brkga(data, decode!, MAXIMIZE, seed, m, brkga_params, perform_evolution)

    # To inject the initial tour, we need to create chromosome representing that
    # solution. First, we create a set of keys to be used in the chromosome.
    # Then, we visit each required and profitable arc in the routes and assign to they keys.
    visited_blocks_weights::Dict{Vi, Real} = Dict{Vi, Real}(visited_blocks[i] => 1 - i * 1e-3 for i::Int in 1:length(visited_blocks))
    initial_chromosome::Vector{Float64}    = [block in keys(visited_blocks_weights) ? visited_blocks_weights[block] : 0.0 for block::Vi in B]

    # Inject the warm start solution in the initial population.
    set_initial_population!(brkga_data, [initial_chromosome])

    # NOTE: don't forget to initialize the algorithm.
    verbose && println("\n[$(Dates.Time(Dates.now()))] Initializing BRKGA data...")
    initialize!(brkga_data)

    ########################################
    # Warm up the script/code
    ########################################

    # To make sure we are timing the runs correctly, we run some warmup
    # iterations with bogus data. Warmup is always recommended for script
    # languages. Here, we call the most used methods.
    verbose && println("\n[$(Dates.Time(Dates.now()))] Warming up...")

    bogus_data = deepcopy(brkga_data)
    evolve!(bogus_data, 2)
    path_relink!(bogus_data, brkga_params.pr_type, brkga_params.pr_selection, (x, y) -> 1.0, (x, y) -> true, 0, 0.5, 1, 10.0, 1.0)
    best_cost::Float64 = get_best_fitness(brkga_data)
    best_chromosome::Vector{Float64} = get_best_chromosome(brkga_data)
    bogus_data = nothing

    ########################################
    # Evolving
    ########################################

    verbose && println("\n[$(Dates.Time(Dates.now()))] Evolving...")
    verbose && println("* Iteration | Cost | CurrentTime")

    iteration::Int = 0
    last_update_time = 0.0
    last_update_iteration = 0
    large_offset = 0
    path_relink_time = 0.0
    num_path_relink_calls::Int = 0
    num_homogenities::Int = 0
    num_best_improvements::Int = 0
    num_elite_improvements::Int = 0
    run::Bool = true
    start_time = Base.time()

    # Main optimization loop. We evolve one generation at time,
    # keeping track of all changes during such process.
    while run
        iteration += 1

        # Evolves one iteration.
        evolve!(brkga_data)

        # Checks the current results and holds the best.
        fitness = get_best_fitness(brkga_data)
        if fitness > best_cost
            last_update_time = Base.time() - start_time
            update_offset = iteration - last_update_iteration

            if large_offset < update_offset
                large_offset = update_offset
            end

            last_update_iteration = iteration

            best_cost = fitness
            best_chromosome = get_best_chromosome(brkga_data)

            verbose && @printf("* %d | %.0f | %.2f \n", iteration, best_cost, last_update_time)

        end

        iter_without_improvement = iteration - last_update_iteration

        # Here, we call the path relink when the algorithm gets stuck for
        # `exchange_interval` iterations. Obviously, we can use many other ways
        # of hybridization.
        if control_params.exchange_interval > 0 &&
            iter_without_improvement > 0 &&
            (iter_without_improvement % control_params.exchange_interval == 0)

            verbose && println("Performing path relink at $iteration...")
            num_path_relink_calls += 1

            pr_now = Base.time()
            result = path_relink!(
                                  brkga_data,
                                  brkga_params.pr_type,
                                  brkga_params.pr_selection,
                                  kendall_tau_distance,
                                  affect_solution_kendall_tau,
                                  brkga_params.pr_number_pairs,
                                  brkga_params.pr_minimum_distance,
                                  1, # block_size doesn't matter for permutation.
                                  maximum_time - (Base.time() - start_time),
                                  brkga_params.pr_percentage
                                 )

            pr_time = Base.time() - pr_now
            path_relink_time += pr_time

            if result == TOO_HOMOGENEOUS
                num_homogenities += 1
                verbose && println("- Populations are too too homogeneous | " *
                                         "Elapsed time: $(@sprintf("%.2f", pr_time))")

            elseif result == NO_IMPROVEMENT
                verbose && println("- No improvement found | " *
                                         "Elapsed time: $(@sprintf("%.2f", pr_time))")

            elseif result == ELITE_IMPROVEMENT
                num_elite_improvements += 1
                verbose && println("- Improvement on the elite set but " *
                                         "not in the best individual | " *
                                         "Elapsed time: $(@sprintf("%.2f", pr_time))")

            elseif result == BEST_IMPROVEMENT
                num_best_improvements += 1
                fitness = get_best_fitness(brkga_data)
                verbose && println("- Best individual improvement: $fitness | " *
                                         "Elapsed time: $(@sprintf("%.2f", pr_time))")
                if fitness < best_cost
                    last_update_time = Base.time() - start_time
                    update_offset = iteration - last_update_iteration

                    if large_offset < update_offset 
                        large_offset = update_offset
                    end

                    last_update_iteration = iteration
                    best_cost = fitness
                    best_chromosome = get_best_chromosome(brkga_data)

                    verbose && @printf("* %d | %.0f | %.2f \n", iteration, best_cost,
                                       last_update_time)
                end
            end
        end

        # Check stop criteria.
        run = !(
                (Base.time() - start_time > maximum_time) ||
                (stop_rule == GENERATIONS && Float64(iteration) == stop_argument) ||
                (stop_rule == IMPROVEMENT &&
                 Float64(iter_without_improvement) >= stop_argument) ||
                (stop_rule == TARGET && best_cost <= stop_argument)
               )
    end
    total_elapsed_time = Base.time() - start_time
    total_num_iterations = iteration

    if verbose
        println("[$(Dates.Time(Dates.now()))] End of optimization\n")
        println("Total number of iterations: $total_num_iterations")
        println("Last update iteration: $last_update_iteration")
        @printf("Total optimization time: %.2f\n", total_elapsed_time)
        @printf("Last update time: %.2f\n", last_update_time)
        println("Large number of iterations between improvements: $large_offset")

        @printf("Total path relink time: %.2f\n", path_relink_time)
        println("Total path relink calls: $num_path_relink_calls")
        println("Number of homogenities: $num_homogenities")
        println("Improvements in the elite set: $num_elite_improvements")
        println("Best individual improvements: $num_best_improvements")
    end

    ########################################
    # Extracting the solution
    ########################################

    # get blocks indexes
    idxs_blocks = getIdxsBlocks(best_chromosome, data)

    # get dijkstra time matrix
    tour_time, consumed_times = dijkstra(data, idxs_blocks)

    fl = findlast(idx::Int -> any(i::Int -> consumed_times[(idx, i)] <= data.T, data.B[idx]), idxs_blocks)
    if fl === nothing || tour_time > data.T
        tour = Vi([data.depot])
        blocks = VVi()
    else
        idxs_feas::Vi = idxs_blocks[begin:fl]
        tour = getDijkstraRoute(data, consumed_times, idxs_feas)
        blocks = VVi([B[idx_block] for idx_block in idxs_feas])
    end

    if fl !== nothing
        println("Last index ", fl)
    else
        println("Last index ", nothing, " (no feasible block prefix under decode / budget T)")
    end

    # log
    info["cost"], info["solverTime"] = string(best_cost), string(total_elapsed_time)

    println("BRKGA cost:", best_cost)
    println("Blocks cost:", sum(data.profits[block] for block in blocks))
    println("Calculated cost:", info["cost"])

    return SBRPSolution(tour, blocks), info

end
