using BrkgaMpIpr

mutable struct SBRPData <: AbstractInstance
    D::InputDigraph
    depot::Int
    B::Vector{Vi}
    T::Float64
    profits::Dict{Vi,Float64}
end

mutable struct SBRPSolution # SBRP Solution
    tour::Vi                 # route
    B::VVi                   # serviced blocks
end

# 40 km / 60 min = 40 km/h
const NORMAL_SPEED::Float64 = (40.0 * 1e3) / 60.0

# arc time in 40 km/h
const time(data::SBRPData, (i, j)::Arc)::Float64 = i == j ? 0.0 : data.D.distance[Arc(i, j)] / NORMAL_SPEED

# arc distance
const distance(data::SBRPData, (i, j)::Arc)::Float64 = i == j ? 0.0 : data.D.distance[Arc(i, j)]

# block distance in meters (length ≤ 1: no perimeter arcs; avoids empty sum and missing closing arc)
const blockDistance(data::SBRPData, block::Vi)::Float64 =
    length(block) <= 1 ? 0.0 :
    (
        sum(
            a::Arc -> data.D.distance[a],
            Arcs(collect(zip(block[begin:end-1], block[begin+1:end])));
            init=0.0,
        ) + data.D.distance[Arc(last(block), first(block))]
    )

const blockDistance(data::SBRPData, distances::ArcCostMap, block::Vi)::Float64 =
    length(block) <= 1 ? 0.0 :
    (
        sum(
            a::Arc -> distances[a],
            Arcs(collect(zip(block[begin:end-1], block[begin+1:end])));
            init=0.0,
        ) + distances[Arc(last(block), first(block))]
    )

# block service time made in 10 km/h
const blockTime(data::SBRPData, block::Vi)::Float64 = 4 * blockDistance(data, block) / NORMAL_SPEED

const blockTime(data::SBRPData, distances::ArcCostMap, block::Vi)::Float64 = 4 * blockDistance(data, distances, block) / NORMAL_SPEED

# distance of a tour
const tourDistance(data::SBRPData, tour::Vi)::Float64 = sum(a::Arc -> data.D.distance[a], Arcs(collect(zip(tour[begin:end-1], tour[begin+1:end]))); init=0.0)

const tourDistance(distances::ArcCostMap, tour::Vi)::Float64 = sum(a::Arc -> distances[a], Arcs(collect(zip(tour[begin:end-1], tour[begin+1:end]))); init=0.0)

# tour time considering the time to spraying all the blocks
const tourTime(data::SBRPData, solution::SBRPSolution)::Float64 = (tourDistance(data, solution.tour) / NORMAL_SPEED) + sum(block::Vi -> blockTime(data, block), solution.B; init=0.0)

const tourTime(data::SBRPData, distances::ArcCostMap, solution::SBRPSolution)::Float64 = (tourDistance(distances, solution.tour) / NORMAL_SPEED) + sum(block::Vi -> blockTime(data, distances, block), solution.B; init=0.0)

#=
Get the number of nodes that appear more than once in a tour
input:
- tour::Vi is the tour
output:
- Number of repeated nodes
=#
function getNumRepeatedNodes(tour::Vi)::Int
    counts = Dict{Int,Int}()
    for node in tour
        counts[node] = get(counts, node, 0) + 1
    end

    repeated_nodes = 0
    for count in values(counts)
        if count > 1
            repeated_nodes += 1
        end
    end

    return repeated_nodes
end

#=
Get the number of arcs that appear more than once in a tour
input:
- tour::Vi is the tour
output:
- Number of repeated arcs
=#
function getNumRepeatedArcs(tour::Vi)::Int
    arc_counts = Dict{Arc,Int}()
    for i in 1:(length(tour)-1)
        arc = Arc(tour[i], tour[i+1])
        arc_counts[arc] = get(arc_counts, arc, 0) + 1
    end

    repeated_arcs = 0
    for count in values(arc_counts)
        if count > 1
            repeated_arcs += 1
        end
    end

    return repeated_arcs
end

#=
Count the occurrences of nodes in a tour that belong to a given set of blocks.
input:
- tour::Vi is the tour.
- serviced_blocks::VVi is the list of blocks.
output:
- The total count of nodes in the tour that are part of the serviced blocks.
=#
function count_nodes_in_serviced_blocks(tour::Vi, serviced_blocks::VVi)::Int
    if isempty(serviced_blocks) || isempty(tour)
        return 0
    end

    serviced_block_nodes = Set(node for block in serviced_blocks for node in block)

    count = 0
    for node in tour
        if node in serviced_block_nodes
            count += 1
        end
    end
    return count
end


# get nodes belonging to some block
const getBlocksNodes(data::SBRPData)::Si = Si(reduce(vcat, data.B))

#=
Drop blocks whose profit is exactly zero from `data.B` and `data.profits`.
Blocks with profit < 0 or > 0 are kept (only `== 0` is removed).
Vertices that belonged only to removed blocks are no longer in any block and
therefore disappear from the compact complete digraph when `createCompleteDigraph`
runs; the street digraph `data.D` is unchanged so shortest-path closure still sees
all Steiner nodes.

Throws `ArgumentError` if no blocks remain or no block has profit > 0 after filtering.
=#
function dropZeroProfitBlocks!(data::SBRPData)::Nothing
    kept::VVi = VVi()
    for block::Vi in data.B
        p::Float64 = data.profits[block]
        if p == 0.0
            delete!(data.profits, block)
        else
            push!(kept, block)
        end
    end
    empty!(data.B)
    append!(data.B, kept)
    if isempty(data.B) || !any(b::Vi -> data.profits[b] > 0, data.B)
        throw(ArgumentError("dropZeroProfitBlocks!: no blocks with profit > 0 remain after removing zero-profit blocks"))
    end
    return nothing
end

# add dummy arcs connected to the nodes belonging to some block
function addDummyArcs(data::SBRPData)

    @debug "Adding dummy arcs"

    for i::Int in getBlocksNodes(data)

        @debug ">> For node $i"

        data.D.distance[Arc(data.depot, i)] = data.D.distance[Arc(i, data.depot)] = 0.0
    end
end

#=
Compute one source row for metric closure: distances and vertex paths from `i` to each
reachable `j ∈ Vb` on the street digraph. Returns fresh `Dict`s for thread-safe merge.
=#
function _completeDigraphRowFromSource(
    i::Int,
    data::SBRPData,
    adjList::Dict{Int,Vi},
    V::Vi,
    Vb::Si,
)::Tuple{Dict{Arc,Float64},Dict{Arc,Vi}}

    distances_local::Dict{Int,Float64} = Dict{Int,Float64}(map(v::Int -> v => typemax(Float64), V))
    pred::Dict{Int,Int} = Dict{Int,Int}(map(v::Int -> v => v, V))
    q::Vi = Vi([i])

    distances_local[i] = 0.0

    while !isempty(q)

        curr::Int = popfirst!(q)

        for next::Int in adjList[curr]

            if !in(next, [data.depot, i]) && distances_local[next] > distances_local[curr] + data.D.distance[Arc(curr, next)]

                distances_local[next] = distances_local[curr] + data.D.distance[Arc(curr, next)]
                pred[next] = curr
                push!(q, next)

            end

        end

    end

    row_dist::Dict{Arc,Float64} = Dict{Arc,Float64}()
    row_paths::Dict{Arc,Vi} = Dict{Arc,Vi}()

    for j::Int in Vb

        pred[j] == j && continue

        row_dist[Arc(i, j)] = distances_local[j]
        path_j::Vi = Vi()
        curr_j::Int = j

        while pred[curr_j] != curr_j

            pushfirst!(path_j, curr_j)
            curr_j = pred[curr_j]

        end

        pushfirst!(path_j, i)
        row_paths[Arc(i, j)] = path_j

    end

    return row_dist, row_paths
end

#=
Create a SBRP instance in which the digraph is complete and metric
input:
- data::SBRPData is a SBRP instance
output:
- data´::SBRPData is the new SBRP instance
- paths´::Dict{Pair{Int, Int}, Vector{Int}} is the relation of paths
=#
function createCompleteDigraph(data::SBRPData)::Pair{SBRPData,Dict{Arc,Vi}}

    @debug "Creating complete digraph"

    # setup
    A::Arcs = data.D.A
    arcs_set::ArcsSet = ArcsSet(data.D.A)
    B::VVi = data.B
    Vb::Si = Si([i for b in B for i in b])
    depot::Int = data.depot

    data′::SBRPData = SBRPData(
        InputDigraph(
            Dict{Int,Vertex}(vcat(
                depot => data.D.V[depot],
                map(i::Int -> i => data.D.V[i], collect(Vb))
            )),
            Arcs(),
            ArcCostMap(map(
                a::Arc -> a => 0.0,
                filter((i, j)::Arc -> i != j, χ(Vb))
            ))
        ),
        depot,
        B,
        data.T, # 2 hours
        data.profits
    )
    V::Vi = Vi(collect(keys(data.D.V)))
    paths::Dict{Arc,Vi} = Dict{Arc,Vi}(map((i, j)::Arc -> Arc(i, j) => Vi(), filter((i, j)::Pair{Int,Int} -> i != j, χ(Vb))))

    # adjList
    @debug "Creating adjacency list"

    adjList::Dict{Int,Vi} = Dict{Int,Vi}(map(i::Int -> i => Vi(), V))

    for (i::Int, j::Int) in A
        push!(adjList[i], j)
    end

    # subrotine to find all shortest paths with root at node i
    function findPaths(i::Int)

        @debug ">> Finding paths from node $i"

        # bfs
        distances::Dict{Int,Float64} = Dict{Int,Float64}(map(i::Int -> i => typemax(Float64), V))
        pred::Dict{Int,Int} = Dict{Int,Int}(map(i::Int -> i => i, V))
        q::Vi = Vi([i])

        distances[i] = 0.0

        while !isempty(q)

            curr::Int = popfirst!(q)

            for next::Int in adjList[curr]

                if !in(next, [data.depot, i]) && distances[next] > distances[curr] + data.D.distance[Arc(curr, next)]

                    distances[next] = distances[curr] + data.D.distance[Arc(curr, next)]
                    pred[next] = curr
                    push!(q, next)

                end

            end

        end

        # add cost and distance
        for j::Int in Vb

            pred[j] == j && continue

            data′.D.distance[Arc(i, j)] = distances[j]
            paths[Arc(i, j)] = Vi()
            curr::Int = j

            while pred[curr] != curr

                pushfirst!(paths[Arc(i, j)], curr)
                curr = pred[curr]

            end

            pushfirst!(paths[Arc(i, j)], i)

        end
    end

    # get paths in parallel
    @debug "Finding paths in parallel"

    path_tasks::Vector{Task} = Vector{Task}()
    for i::Int in Vb
        t = @task findPaths(i)
        schedule(t)
        push!(path_tasks, t)
    end

    # run tasks
    for t::Task in path_tasks
        wait(t)
    end

    # check feasibility of the generated paths
    for (_, path::Vi) in paths

        path_arcs::Vector{Tuple{Int,Int}} = collect(zip(path[begin:end-1], path[begin+1:end]))

        idx::Union{Nothing,Int} = findfirst((i, j)::Tuple{Int,Int} -> !in(Arc(i, j), arcs_set), path_arcs)

        if idx != nothing
            throw(InvalidStateException("In path $path, the arc $((path[idx], path[idx + 1])) does not exist", :closed))
        end

    end

    # dummy arcs
    @debug "Adding dummy arc weights"

    for i::Int in Vb
        data′.D.distance[Arc(i, depot)] = data′.D.distance[Arc(depot, i)] = 0.0
    end

    # add arcs
    data′.D.A = collect(keys(data′.D.distance))

    # remove beginnings and endings of each path
    for (a::Arc, path::Vi) in paths

        popfirst!(path)
        pop!(path)

    end

    # check feasibility
    checkSBRPCompleteFeasibility(data′, Vb)


    return Pair{SBRPData,Dict{Arc,Vi}}(data′, paths)
end

#=
Ordered vertex list for one disease block from its street arcs (mutates `arcs`).
=#
function _carlosBlockVerticesFromArcs!(arcs::Arcs)::Vi
    block::Vi = Vi()

    curr::Int, next::Int = first(arcs)

    push!(block, curr)
    deleteat!(arcs, 1)

    while true
        idx = findfirst((i, j)::Pair{Int,Int} -> next == i, arcs)
        idx == nothing && break

        curr, next = arcs[idx]

        push!(block, curr)
        deleteat!(arcs, idx)
    end

    return block
end

#=
Shared post-read pipeline for Carlos-format SBRP (reverse arcs, depot, metric closure).
`data` must already have `D.V`, `D.distance`, `B`, and `profits` from the file (before reverse arcs).
=#
function _finalizeCarlosSbrpData(
    data::SBRPData,
    app::Dict{String,Any},
)::Tuple{SBRPData,Union{Dict{Arc,Vi},Nothing},Union{ArcCostMap,Nothing}}

    # Undirected streets: add reverse arcs when missing (sparse Path-CBRP / C++ Carlos parity).
    for a::Arc in collect(keys(data.D.distance))
        rev_a::Arc = Arc(last(a), first(a))
        if !haskey(data.D.distance, rev_a)
            data.D.distance[rev_a] = data.D.distance[a]
        end
    end

    # update arcs ids
    data.D.A = collect(keys(data.D.distance))

    # add dummy depot
    @debug "Finding dummy depot"

    depot::Int = data.depot = maximum(collect(keys(data.D.V))) + 1
    data.D.V[depot] = Vertex(depot, -1, -1)

    # dummy weights
    addDummyArcs(data)

    data.D.A = collect(keys(data.D.distance))

    if get(app, "drop-zero-profit-blocks", false)
        dropZeroProfitBlocks!(data)
    end

    # check feasibility
    @debug "Checking instance connectivity"

    Vb::Si = getBlocksNodes(data)
    distances::ArcCostMap = calculateShortestPaths(data)

    not_connected::Union{Nothing,Arc} = nothing
    for p in χ(Vb)
        a_key::Arc = Arc(p.first, p.second)
        if !haskey(distances, a_key) || distances[a_key] >= ∞ / 2
            not_connected = a_key
            break
        end
    end
    if not_connected !== nothing
        throw(ArgumentError("The SBRP instance it is not connected (first missing pair: $(not_connected))"))
    end

    if get(app, "no-cbrp-metric-closure", false)
        @debug "Skipping metric closure (sparse Carlos digraph)"
        return data, nothing, nothing
    end

    # Same object Matheus `readSBRPData` passes to the complete digraph IP: metric-complete digraph
    # on depot ∪ block nodes only (not the street hull `data‴` with full `D.V` and sparse `D.A`).
    data′::SBRPData, paths′::Dict{Arc,Vi} = createCompleteDigraph(data)

    # check feasibility
    checkArcsFeasibility(data, data′)

    street_distances::ArcCostMap = deepcopy(data.D.distance)
    return data′, paths′, street_distances
end

#=
Create a SBRP instance from a Carlo's instance
input:
- app::Dict{String, Any} is the paramenters relation
output:
- `data`: sparse street digraph if `no-cbrp-metric-closure`; otherwise compact metric-complete digraph (`data′`, same as Matheus complete IP input).
- `paths` / `street_distances`: populated iff metric closure ran; used to expand a compact IP tour back to street vertices (`retrieveOriginalDigraphSolution`), same role as `readSBRPData` for Matheus.
=#
function readSBRPDataCarlos(app::Dict{String,Any})::Tuple{SBRPData,Union{Dict{Arc,Vi},Nothing},Union{ArcCostMap,Nothing}}

    @debug "Reading Carlos SBRP instances"

    data::SBRPData = SBRPData(
        InputDigraph(
            Dict{Int,Vertex}(1 => Vertex(-1, -1, -1)),
            Arcs(),
            ArcCostMap()
        ),
        -1,
        VVi(),
        parse(Float64, app["vehicle-time-limit"]),
        ArcCostMap()
    )

    open(app["instance"]) do f::IOStream

        nNodes::Int, nArcs::Int, nBlocks::Int = map(str::String -> parse(Int, str), Vector{String}(split(readline(f), [' ']; limit=3, keepempty=false)))

        # get nodes
        @debug "Reading nodes"

        for i::Int in 1:nNodes

            strs::Vector{String} = split(readline(f), [' ']; keepempty=false)

            id::Int = parse(Int, strs[2])
            x::Float64 = parse(Float64, strs[3])
            y::Float64 = parse(Float64, strs[4])

            data.D.V[id] = Vertex(id, x, y)
        end

        # get arcs
        @debug "Reading arcs"

        block_arcs::Dict{Int,Arcs} = Dict{Int,Arcs}()
        block_profit::Dict{Int,Int} = Dict{Int,Int}()

        for n::Int in 1:nArcs

            strs::Vector{String} = split(readline(f), [' ']; keepempty=false)

            a::Arc = Arc(parse(Int, strs[2]), parse(Int, strs[3]))
            distance::Float64 = parse(Float64, strs[4])
            block_id::Int = parse(Int, strs[5])
            cases::Int = parse(Int, strs[6])

            #rounded_distance = ceil(distance * 10.0) / 10.0
            #rounded_distance = floor(distance)

            data.D.distance[a] = distance

            block_id == -1 && continue

            # new block
            if !haskey(block_profit, block_id)
                block_profit[block_id] = 0
                block_arcs[block_id] = Arcs()
            end

            block_profit[block_id] += cases

            push!(block_arcs[block_id], a)

        end

        # get blocks
        @debug "Reading blocks"

        @debug "Reading blocks profits"
        for (block_id::Int, arcs::Arcs) in block_arcs

            block::Vi = _carlosBlockVerticesFromArcs!(arcs)
            push!(data.B, block)

            data.profits[block] = Float64(block_profit[block_id])
        end

    end

    return _finalizeCarlosSbrpData(data, app)
end

#=
Create a SBRP instance from a Matheus's instance
input:
- app::Dict{String, Any} is the paramenters relation
output:
- data´::SBRPData is the SBRP instance
- paths´::Dict{Arc, Vi} is the relation of paths
- distances::ArcCostMap is the original arc distances relationship
=#
function readSBRPData(app::Dict{String,Any})::Tuple{SBRPData,Dict{Arc,Vi},ArcCostMap}

    @debug "Reading Matheus's instance"

    depot::Int = 1
    data::SBRPData = SBRPData(
        InputDigraph(
            Dict{Int,Vertex}(1 => Vertex(-1, -1, -1)),
            Arcs(),
            ArcCostMap()
        ),
        depot,
        VVi(),
        parse(Float64, app["vehicle-time-limit"]),
        ArcCostMap()
    )
    blocks::Dict{Int,Arcs} = Dict{Int,Arcs}()

    open(app["instance"]) do f::IOStream

        # ignore the first two lines
        readline(f)
        readline(f)

        nNodes::Int = parse(Int, split(readline(f), [' ']; limit=0, keepempty=false)[end])
        nArcs::Int = parse(Int, split(readline(f), [' ']; limit=0, keepempty=false)[end]) + parse(Int, split(readline(f), [' ']; limit=0, keepempty=false)[end])
        nBlocks::Int = parse(Int, split(readline(f), [' ']; limit=0, keepempty=false)[end])

        readline(f)
        readline(f)

        # get arcs
        @debug "Reading arcs"

        for _ in 1:nArcs

            parts::Vector{String} = split(readline(f), ['(', ')', ' ', ',']; limit=0, keepempty=false)

            data.D.distance[Arc(parse(Int, parts[1]), parse(Int, parts[2]))] = floor(Int, parse(Float64, parts[end]))
        end

        # get nodes
        @debug "Reading nodes"

        readline(f)
        for _ in 1:nNodes
            parts::Vector{String} = split(readline(f), [' ']; limit=0, keepempty=false)
            id::Int = parse(Int, parts[1])
            data.D.V[id] = Vertex(id, parse(Float64, parts[2]), parse(Float64, parts[3]))
        end

        # get blocks
        @debug "Reading blocks"

        readline(f)
        for _ in 1:nBlocks
            parts::Vector{String} = split(readline(f), [',', ' ']; limit=0, keepempty=false)
            block::Vi = Vi([parse(Int, part) for part in parts[begin:end-1]])
            push!(data.B, block)

            # define profit
            if app["cluster-size-profits"]
                data.profits[block] = length(block)
            elseif app["unitary-profits"]
                data.profits[block] = 1.0
            else
                data.profits[block] = parse(Float64, parts[end])
            end

        end

    end

    # update arcs ids
    data.D.A = collect(keys(data.D.distance))

    # dummy weights
    addDummyArcs(data)

    if get(app, "drop-zero-profit-blocks", false)
        dropZeroProfitBlocks!(data)
    end

    # check feasibility
    @debug "Checking connectivity"

    Vb::Si = getBlocksNodes(data)

    distances::ArcCostMap = calculateShortestPaths(data)

    if any(a::Arc -> !in(a, keys(distances)), χ(Vb))
        throw(InvalidStateException("The SBRP instance it is not connected"))
    end

    # compact in complete graph
    data′::SBRPData, paths′::Dict{Arc,Vi} = createCompleteDigraph(data)

    # check feasibility
    checkArcsFeasibility(data, data′)

    # consider only shortest paths arcs
    # copy
    data‴ = deepcopy(data)

    # set arcs
    data‴.D.A = collect(ArcsSet(vcat(
        [Arc(i, j) for (a::Arc, path::Vi) in paths′ for (i::Int, j::Int) in zip(path[1:end-1], path[2:end])], # min paths arcs
        [Arc(a.first, path[begin]) for (a::Arc, path::Vi) in paths′ if !isempty(path)], # min paths arcs
        [Arc(path[end], a.second) for (a::Arc, path::Vi) in paths′ if !isempty(path)], # min paths arcs
        [a for (a::Arc, path::Vi) in paths′ if isempty(path)], # min paths arcs (edge case)
        [Arc(i, j) for block::Vi in data.B for (i::Int, j::Int) in zip(block[1:end-1], block[2:end])], # blocks arcs
        map(block::Vi -> Arc(block[end], block[begin]), data.B), # blocks arcs
        map(i::Int -> Arc(depot, i), collect(Vb)), # depot arcs
        map(i::Int -> Arc(i, depot), collect(Vb))  # depot arcs
    )))
    data‴.D.distance = ArcCostMap(map(a::Arc -> a => data.D.distance[a], data‴.D.A))

    # check feasibility
    checkArcsFeasibility(data, data‴)

    # set
    distance = deepcopy(data.D.distance)
    data = data‴

    # return
    return data′, paths′, distance
end

#=
Check arcs match of a SBRP instance: This functions checks if the shortest paths of both instances match
input:
- data::SBRPData is the SBRP instance used for basis
- data´::SBRPData is the SBRP that will be verified
output:
- None
exception:
- Case some irregularity is found
=#
function checkArcsFeasibility(data::SBRPData, data′::SBRPData)

    Vb::Si = getBlocksNodes(data)

    @debug "Calculating shortest paths for the original instance"
    distances::ArcCostMap = calculateShortestPaths(data)

    @debug "Calculating shortest paths for the new instance"
    distances′::ArcCostMap = calculateShortestPaths(data′)

    if any(a::Arc -> !in(a, keys(distances′)), χ(Vb))
        throw(InvalidStateException("The SBRP instance it is not connected", :closed))
    end

    for a::Arc in χ(Vb)
        if abs(distances[a] - distances′[a]) > 1e-3
            throw(InvalidStateException("Original digraph distance of arc $a is $(distances[a]) while the new digraph distance is $(distances′[a])", :closed))
        end
    end
end

#=
Calculate all pair of shortest paths
input:
- data::SBRPData is the SBRP instance
output:
- distances::ArcCostMap is the distances relationship
=#
function calculateShortestPaths(data::SBRPData)::ArcCostMap

    # params
    V::Dict{Int,Vertex} = data.D.V
    A::Arcs = data.D.A
    Vb::Si = getBlocksNodes(data)

    # building adjacency list
    @debug "Building adjacency list"

    adjList::Dict{Int,Vi} = Dict{Int,Vi}(map(i::Int -> i => Vi(), collect(keys(V))))

    for (i::Int, j::Int) in A
        push!(adjList[i], j)
    end

    distances::ArcCostMap = ArcCostMap(map(
        a::Arc -> a => typemax(Float64),
        reduce(vcat, map(i::Int -> map(j::Int -> Arc(i, j), collect(keys(V))), collect(Vb)))
    ))


    # subrotine to find all shortest paths with root at node i
    function findPaths(i::Int)

        @debug ">> Finding paths from node $i"

        # bfs
        q::Vi = Vi([i])

        distances[Arc(i, i)] = 0.0

        while !isempty(q)

            # curr node
            curr::Int = popfirst!(q)

            # next node
            for next::Int in adjList[curr]

                # cost
                b::Float64 = distances[Arc(i, curr)] + data.D.distance[Arc(curr, next)]

                # update distance and push new nodes
                if !in(next, [i, data.depot]) && distances[Arc(i, next)] > b

                    distances[Arc(i, next)] = b
                    push!(q, next)

                end

            end

        end

    end

    # get paths
    @debug "Obtaining shortest paths in parallel"

    path_tasks::Vector{Task} = Vector{Task}()

    for i::Int in Vb
        t::Task = @task findPaths(i)
        schedule(t)
        push!(path_tasks, t)
    end


    for t::Task in path_tasks
        wait(t)
    end

    return distances
end

#=
Check feasibility of a SBRP complete instance
input:
- data::SBRPData is the SBRP instance
- Vb::Set{Int} is the set of nodes belonging to a block
output:
- None
exception:
- Case some irregularity is found
=#
function checkSBRPCompleteFeasibility(data::SBRPData, Vb::Si)

    @debug "Checking feasibility of the SBRP instance on complete digraph"

    if any(a::Arc -> !in(a, data.D.A), filter((i, j)::Arc -> i != j, χ(Vb)))
        throw(InvalidStateException("The Complete SBRP instance is not feasible"))
    end
end


#=
Check feasibility of a SBRP instance
input:
- data::SBRPData is the SBRP instance
output:
- None
exception:
- Case some irregulatiry in data´ is found
=#
function checkSBRPfeasibility(data::SBRPData)

    # params
    V′::Si = getBlocksNodes(data)
    A::Arc = data.D.A
    adjList::Dict{Int,Vi} = Dict{Int,Vi}(map(i::Int -> i => Vi(), V′))
    revAdjList::Dict{Int,Vi} = Dict{Int,Vi}(map(i::Int -> i => Vi(), V′))

    # building adjacency list
    @debug "Building adjacency list"

    for (i::Int, j::Int) in A

        push!(adjList[i], j)
        push!(revAdjList[j], i)

    end

    # Kosaraju’s algorithm
    @debug "Kosaraju’s algorithm"

    # outgoing arcs
    @debug "Outgoing arc"

    # bfs
    v::Int = first(V′)
    mandatory_visited_vertices::Si = Si([v])
    visited_vertices::Si = Si([v])
    q::Vi = Vi([v])

    while !isempty(q)

        curr::Int = popfirst!(q)

        for j::Int in adjList[curr]

            if !in(j, visited_vertices)

                push!(q, j)
                push!(visited_vertices, j)

                if j::Int in V′

                    push!(mandatory_visited_vertices, j)

                    if length(mandatory_visited_vertices) == length(V′)

                        empty(q)
                        break

                    end
                end
            end
        end
    end

    if length(mandatory_visited_vertices) < length(V′)
        for u::Int in V′
            if !in(u, mandatory_visited_vertices)
                throw(InvalidStateException("There is not path between $v and $u"))
            end
        end
    end

    # ingoing arcs
    @debug "Ingoing arc"

    # bfs
    mandatory_visited_vertices = Si([v])
    visited_vertices = Si([v])
    q = [v]
    while !isempty(q)

        curr::Int = popfirst!(q)

        for i::Int in revAdjList[curr]

            if !in(i, visited_vertices)

                push!(q, i)
                push!(visited_vertices, i)

                if i::Int in V′

                    push!(mandatory_visited_vertices, i)

                    if length(mandatory_visited_vertices) == length(V′)
                        empty(q)
                        break
                    end
                end
            end
        end
    end

    if length(mandatory_visited_vertices) < length(V′)
        for u::Int in V′
            if !in(u, mandatory_visited_vertices)
                throw(InvalidStateException("There is not path between $v and $u"))
            end
        end
    end
end
