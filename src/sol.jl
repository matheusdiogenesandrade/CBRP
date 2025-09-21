#= 
Write solution in a file
input:
- file_path::String is the file path the solution file will be saved
- tour::Vector{Vertex} is the route to be written in a file
- data::SBRPData is the SBRP instance in which the tour was built on
- B::Vector{Vector{Int}} is the set of blocks serviced by the route ´route´
output:
- None
=#
function writeSolution(file_path::String, data::SBRPData, solution::SBRPSolution)

    @debug "Writing solution in a file"

    # get solution params
    tour::Vi = solution.tour

    serviced_blocks::VVi = solution.B

    # write route
    tour_no_depot::Vi = filter(i::Int -> i != data.depot, tour)

    sol_string::String = join(tour_no_depot, ", ") * "\n"

    for block::Vi in serviced_blocks
        sol_string *= join(block, ", ") * "\n"
    end

    write(open(file_path * ".sol", "w"), sol_string)

    # write GPS
#    writeGPX(file_path * ".gpx", map(i::Int -> data.D.V[i], tour_no_depot))
end

#= 
Check SBRP solution feasibility
input:
- data::Vector{Vertex} is the route to be written in a file
output:
- None
=#
function checkSBRPSolution(data::SBRPData, solution::SBRPSolution)

    @debug "Checking the feasibility of a SBRP solution"

    # instance parameters
    A::ArcsSet = ArcsSet(data.D.A)

    # solution parameters
    tour::Vi = solution.tour
    visited_nodes::Si = Set{Int}(tour)
    serviced_blocks::VVi = solution.B

    # check arcs
    for (i::Int, j::Int) in zip(tour[begin:end - 1], tour[begin + 1:end])

        if !in(Arc(i, j), A) 
            throw(ErrorException("Arc $(Arc(i, j)) does not exists"))
        end

    end

    # check blocks
    for block::Vi in serviced_blocks

        if all(i::Int -> !in(i, visited_nodes), block)
            throw(ErrorException("Block $block was not served"))
        end

    end

end

#= 
Read a SBRPSolution from a file. If the file contains only a tour and no blocks,
it runs a knapsack algorithm to select a profitable set of blocks along the tour.
input:
- data::SBRPData is the SBRP instance
- file_path::String is the path to the solution file
- distances::Dict{Arc, Float64} contains all-pairs shortest path distances.
output:
- SBRPSolution object
=#
function readSBRPSolution(data::SBRPData, file_path::String, distances::Dict{Arc, Float64})::SBRPSolution
    
    lines = readlines(file_path)
    
    # Parse the tour (first line)
    tour_strings = split(lines[1], ",")
    tour_nodes = [parse(Int, strip(s)) for s in tour_strings if !isempty(strip(s))]
    tour = tour_nodes
    
    # Parse the blocks (remaining lines)
    serviced_blocks = VVi()
    for line in lines[2:end]
        if !isempty(strip(line))
            block_strings = split(line, ",")
            block_nodes = [parse(Int, strip(s)) for s in block_strings if !isempty(strip(s))]
            push!(serviced_blocks, block_nodes)
        end
    end

    # If no blocks were read from the file, run the knapsack algorithm
    if isempty(serviced_blocks) && !isempty(tour)
        @debug "No blocks in solution file, running knapsack to select blocks."

        # 1. Calculate remaining time (Knapsack Capacity)
        # Create a temporary solution object to calculate the tour's time
        temp_solution = SBRPSolution(tour, VVi())
        tour_travel_time = tourTime(data, distances, temp_solution)
        capacity = data.T - tour_travel_time

        if capacity > 0
            # 2. Identify candidate blocks (items)
            tour_nodes_set = Set(tour)
            candidate_blocks = [block for block in data.B if any(node -> node in tour_nodes_set, block)]

            if !isempty(candidate_blocks)
                # 3. Get weights and values
                weights = [blockTime(data, block) for block in candidate_blocks]
                values  = [data.profits[block] for block in candidate_blocks]

                # 4. Solve the knapsack problem
                serviced_blocks = knapsack_solver(capacity, weights, values, candidate_blocks)
                @debug "Knapsack selected $(length(serviced_blocks)) blocks."
            end
        end
    end
    
    return SBRPSolution(tour, serviced_blocks)
end


#= 
Calculate the distance of a tour segment.
input:
- distances::Dict{Arc, Float64} is the dictionary mapping arcs to distances.
- path::Vi is a sequence of nodes.
output:
- The total distance of the path.
=#
function get_path_distance(distances::Dict{Arc, Float64}, path::Vi)::Float64
    dist = 0.0
    if length(path) < 2
        return dist
    end

    for i in 1:(length(path) - 1)
        u, v = path[i], path[i+1]
        dist += distances[Arc(u, v)]
    end
    return dist
end

#= 
Calculate the average detour index for a given solution.
This measures the average inefficiency of travel between service blocks.
An index of 0.0 means all travel between blocks is via a shortest path.
An index of 0.1 means travel is 10% longer than the shortest path on average.
input:
- original_tour::Vi is the fully expanded sequence of nodes.
- solution::SBRPSolution contains the list of visited blocks.
- distances::Dict{Arc, Float64} contains all-pairs shortest path distances.
output:
- The average detour index as a Float64.
=#
function calculate_average_detour_index(
    original_tour::Vi,
    solution::SBRPSolution,
    distances::Dict{Arc, Float64}
)::Float64

    # 1. Create a Set of all nodes in visited blocks for O(1) lookup.
    visited_block_nodes = Set(node for block in solution.B for node in block)

    # 2. Single pass to find inter-block segments and calculate detour.
    total_detour_index = 0.0
    num_segments = 0
    
    current_segment = Vi()
    is_traveling_between_blocks = false

    if isempty(original_tour)
        return 0.0
    end

    # Check if the tour starts outside a block
    if !(original_tour[1] in visited_block_nodes)
        is_traveling_between_blocks = true
    end

    for i in 1:(length(original_tour) - 1)
        u = original_tour[i]
        v = original_tour[i+1]

        u_in_block = u in visited_block_nodes
        v_in_block = v in visited_block_nodes

        # State: Continuing travel between blocks
        if is_traveling_between_blocks && !v_in_block
            push!(current_segment, u)
        
        # State transition: Exiting a block
        elseif u_in_block && !v_in_block
            is_traveling_between_blocks = true
            current_segment = [u] # Start a new segment
        
        # State transition: Entering a block
        elseif is_traveling_between_blocks && v_in_block
            is_traveling_between_blocks = false
            push!(current_segment, u)
            push!(current_segment, v)

            # --- Process the completed segment ---
            if length(current_segment) >= 2
                segment_start_node = current_segment[1]
                segment_end_node = current_segment[end]

                # 1. Calculate actual distance of the segment
                actual_dist = get_path_distance(distances, current_segment)

                # 2. Look up optimal distance (O(1))
                optimal_dist = distances[Arc(segment_start_node, segment_end_node)]

                # 3. Accumulate detour index
                if optimal_dist > 0
                    total_detour_index += (actual_dist / optimal_dist) - 1
                    num_segments += 1
                end
            end
            current_segment = Vi() # Reset
        end
    end

    if num_segments == 0
        return 0.0 # No inter-block segments found
    end

    return total_detour_index / num_segments
end

#= 
Solves the 0/1 knapsack problem using dynamic programming.
input:
- capacity::Float64 is the maximum weight the knapsack can hold.
- weights::Vector{Float64} are the weights of the items.
- values::Vector{Float64} are the values of the items.
- items::VVi is the actual list of items (blocks) to choose from.
output:
- A VVi containing the items selected for the knapsack.
=#
function knapsack_solver(capacity::Float64, weights::Vector{Float64}, values::Vector{Float64}, items::VVi)::VVi
    n = length(items)
    # Note: This DP approach assumes integer weights. We floor the float times to use it.
    # This is a simplification; a more complex knapsack algorithm would be needed for arbitrary float weights.
    int_capacity = floor(Int, capacity)
    int_weights = map(w -> floor(Int, w), weights)

    if int_capacity <= 0 || n == 0
        return VVi()
    end

    # dp[i][w] stores the max value using items 1..i-1 with capacity w-1
    dp = zeros(Float64, n + 1, int_capacity + 1)

    for i in 1:n
        for w in 1:int_capacity
            if int_weights[i] <= w
                dp[i+1, w+1] = max(dp[i, w+1], values[i] + dp[i, w+1 - int_weights[i]])
            else
                dp[i+1, w+1] = dp[i, w+1]
            end
        end
    end

    # Backtrack to find the selected items
    selected_items = VVi()
    res = dp[n+1, int_capacity+1]
    w = int_capacity
    for i in n:-1:1
        if res <= 0 break end
        # if the result comes from including the current item
        if res != dp[i, w+1]
            push!(selected_items, items[i])
            res -= values[i]
            w -= int_weights[i]
        end
    end
    
    return selected_items
end
