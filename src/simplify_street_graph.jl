#=
OSMnx-style topological street digraph simplification (block-safe).

Mirrors `osmnx.simplification._is_endpoint` / `_build_path` / `simplify_graph`
(Boeing 2025), without geometry attributes. Lengths of contracted corridors
are summed. Nodes in `keep` (depot + all block members) are always endpoints.
=#

"""
    streetDigraphNeighbors(A) -> (succ, pred, neighbors)

Adjacency maps for a simple digraph given as an arc list.
"""
function streetDigraphNeighbors(A::Arcs)::Tuple{Dict{Int,Vi},Dict{Int,Vi},Dict{Int,Si}}
    succ::Dict{Int,Vi} = Dict{Int,Vi}()
    pred::Dict{Int,Vi} = Dict{Int,Vi}()
    neighbors::Dict{Int,Si} = Dict{Int,Si}()
    for (u::Int, v::Int) in A
        push!(get!(Vi, succ, u), v)
        push!(get!(Vi, pred, v), u)
        push!(get!(Si, neighbors, u), v)
        push!(get!(Si, neighbors, v), u)
    end
    return succ, pred, neighbors
end

"""
    isStreetEndpoint(node, succ, pred, neighbors, keep) -> Bool

OSMnx `_is_endpoint` rules 1–3, plus forced keep-set (rule 4 analogue).
`degree` is in-degree + out-degree (NetworkX MultiDiGraph.degree).
"""
function isStreetEndpoint(
    node::Int,
    succ::Dict{Int,Vi},
    pred::Dict{Int,Vi},
    neighbors::Dict{Int,Si},
    keep::Si,
)::Bool
    node in keep && return true

    nbrs::Si = get(neighbors, node, Si())
    n::Int = length(nbrs)
    in_d::Int = length(get(pred, node, Vi()))
    out_d::Int = length(get(succ, node, Vi()))
    d::Int = in_d + out_d

    # RULE 1: self-loop
    node in nbrs && return true

    # RULE 2: all inbound or all outbound
    (out_d == 0 || in_d == 0) && return true

    # RULE 3: not (exactly 2 neighbors and degree 2 or 4)
    if !((n == 2) && (d == 2 || d == 4))
        return true
    end

    return false
end

"""
    buildStreetSimplifyPath(endpoint, endpoint_successor, endpoints, succ) -> Vi

OSMnx `_build_path`: walk from `endpoint` through `endpoint_successor` until
the next endpoint. Path endpoints are true endpoints; interior nodes are
interstitial.
"""
function buildStreetSimplifyPath(
    endpoint::Int,
    endpoint_successor::Int,
    endpoints::Si,
    succ::Dict{Int,Vi},
)::Vi
    path::Vi = Vi([endpoint, endpoint_successor])
    path_set::Si = Si(path)

    for this_successor::Int in get(succ, endpoint_successor, Vi())
        successor::Int = this_successor
        if successor in path_set
            continue
        end
        push!(path, successor)
        push!(path_set, successor)
        while !(successor in endpoints)
            successors::Vi = Vi([n for n in get(succ, successor, Vi()) if !(n in path_set)])
            if length(successors) == 1
                successor = successors[1]
                push!(path, successor)
                push!(path_set, successor)
            elseif length(successors) == 0
                if endpoint in get(succ, successor, Vi())
                    push!(path, endpoint)
                    return path
                end
                @warn "Unexpected simplify pattern handled near $successor"
                return path
            else
                throw(ArgumentError("Impossible simplify pattern failed near $successor"))
            end
        end
        return path
    end
    return path
end

"""
    collectStreetSimplifyPaths(succ, pred, neighbors, keep) -> (endpoints, paths)

Identify endpoints and yield every interstitial path to contract.
"""
function collectStreetSimplifyPaths(
    succ::Dict{Int,Vi},
    pred::Dict{Int,Vi},
    neighbors::Dict{Int,Si},
    nodes::Si,
    keep::Si,
)::Tuple{Si,Vector{Vi}}
    endpoints::Si = Si([
        n for n in nodes if isStreetEndpoint(n, succ, pred, neighbors, keep)
    ])
    paths::Vector{Vi} = Vi[]
    for endpoint::Int in endpoints
        for successor::Int in get(succ, endpoint, Vi())
            if !(successor in endpoints)
                push!(paths, buildStreetSimplifyPath(endpoint, successor, endpoints, succ))
            end
        end
    end
    return endpoints, paths
end

"""
    pathCorridorLength(path, distance) -> Float64

Sum of arc lengths along consecutive pairs of `path`.
"""
function pathCorridorLength(path::Vi, distance::ArcCostMap)::Float64
    total::Float64 = 0.0
    for i::Int in 1:(length(path) - 1)
        a::Arc = Arc(path[i], path[i + 1])
        haskey(distance, a) || throw(ArgumentError("Missing arc $a on simplify path $path"))
        total += distance[a]
    end
    return total
end

"""
    weaklyConnectedComponents(nodes, neighbors) -> Vector{Si}
"""
function weaklyConnectedComponents(nodes::Si, neighbors::Dict{Int,Si})::Vector{Si}
    remaining::Si = copy(nodes)
    comps::Vector{Si} = Si[]
    while !isempty(remaining)
        start::Int = first(remaining)
        comp::Si = Si()
        stack::Vi = Vi([start])
        while !isempty(stack)
            u::Int = pop!(stack)
            u in comp && continue
            push!(comp, u)
            delete!(remaining, u)
            for v::Int in get(neighbors, u, Si())
                if v in remaining
                    push!(stack, v)
                end
            end
        end
        push!(comps, comp)
    end
    return comps
end

"""
    validateBlockArcsPresent!(data)

Ensure every consecutive block arc (and closing arc) exists in `data.D.distance`.
"""
function validateBlockArcsPresent!(data::SBRPData)
    for block::Vi in data.B
        length(block) <= 1 && continue
        for (i::Int, j::Int) in zip(block[begin:(end - 1)], block[(begin + 1):end])
            haskey(data.D.distance, Arc(i, j)) ||
                throw(ArgumentError("Block arc missing after simplify: $(Arc(i, j)) in $block"))
        end
        haskey(data.D.distance, Arc(last(block), first(block))) || throw(
            ArgumentError(
                "Block closing arc missing after simplify: $(Arc(last(block), first(block)))",
            ),
        )
    end
    return nothing
end

"""
    simplifyStreetDigraph!(data; remove_rings=true) -> NamedTuple

In-place OSMnx-style simplification. Forces depot and all block member nodes
to remain. Returns size stats for logging.
"""
function simplifyStreetDigraph!(
    data::SBRPData;
    remove_rings::Bool=true,
)::NamedTuple{(:nodes_before, :nodes_after, :arcs_before, :arcs_after),NTuple{4,Int}}
    nodes_before::Int = length(data.D.V)
    arcs_before::Int = length(data.D.A)

    keep::Si = Si(getBlocksNodes(data))
    push!(keep, data.depot)

    for v::Int in keep
        haskey(data.D.V, v) || throw(ArgumentError("Keep-set node $v missing from digraph"))
    end

    succ, pred, neighbors = streetDigraphNeighbors(data.D.A)
    nodes::Si = Si(keys(data.D.V))
    endpoints, paths = collectStreetSimplifyPaths(succ, pred, neighbors, nodes, keep)

    nodes_to_remove::Si = Si()
    new_edges::ArcCostMap = ArcCostMap()
    for path::Vi in paths
        length(path) < 2 && continue
        for interstitial::Int in path[2:(end - 1)]
            interstitial in keep && throw(
                ArgumentError("Refusing to remove keep-set node $interstitial from path $path"),
            )
            push!(nodes_to_remove, interstitial)
        end
        u::Int = path[begin]
        v::Int = path[end]
        len::Float64 = pathCorridorLength(path, data.D.distance)
        a::Arc = Arc(u, v)
        if haskey(new_edges, a)
            new_edges[a] = min(new_edges[a], len)
        else
            new_edges[a] = len
        end
    end

    # Keep arcs whose both ends survive; then merge contracted corridors.
    distance′::ArcCostMap = ArcCostMap()
    for (a::Arc, c::Float64) in data.D.distance
        u::Int, v::Int = a
        if !(u in nodes_to_remove) && !(v in nodes_to_remove)
            distance′[a] = c
        end
    end
    for (a::Arc, c::Float64) in new_edges
        if haskey(distance′, a)
            distance′[a] = min(distance′[a], c)
        else
            distance′[a] = c
        end
    end

    # Drop removed vertices
    for v::Int in nodes_to_remove
        delete!(data.D.V, v)
    end

    if remove_rings
        A_tmp::Arcs = collect(keys(distance′))
        succ′, pred′, neighbors′ = streetDigraphNeighbors(A_tmp)
        nodes′::Si = Si(keys(data.D.V))
        for comp::Si in weaklyConnectedComponents(nodes′, neighbors′)
            if !any(n -> isStreetEndpoint(n, succ′, pred′, neighbors′, keep), comp)
                for n::Int in comp
                    n in keep && continue
                    delete!(data.D.V, n)
                    for a::Arc in collect(keys(distance′))
                        if a.first == n || a.second == n
                            delete!(distance′, a)
                        end
                    end
                end
            end
        end
    end

    # Drop any arc incident to a deleted vertex (ring cleanup / safety)
    alive::Si = Si(keys(data.D.V))
    for a::Arc in collect(keys(distance′))
        if !(a.first in alive) || !(a.second in alive)
            delete!(distance′, a)
        end
    end

    data.D.distance = distance′
    data.D.A = collect(keys(distance′))

    for v::Int in keep
        haskey(data.D.V, v) || throw(ArgumentError("Keep-set node $v missing after simplify"))
    end
    validateBlockArcsPresent!(data)

    return (
        nodes_before=nodes_before,
        nodes_after=length(data.D.V),
        arcs_before=arcs_before,
        arcs_after=length(data.D.A),
    )
end
