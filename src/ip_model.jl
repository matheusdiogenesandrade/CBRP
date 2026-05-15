using Printf
using Combinatorics
using Logging

include("SparseMaxFlowMinCut.jl")

"""CPLEX per-solve cap from `--time-limit` (0 → 3600)."""
function _ip_base_cplex_seconds(app_time_limit::Float64)::Float64
    return app_time_limit > 0.0 ? app_time_limit : 3600.0
end

"""Limit CPLEX global thread count (barrier, parallel MIP, etc.)."""
function _set_cplex_threads!(model::Model, n::Int=1)
    set_optimizer_attribute(model, "CPXPARAM_Threads", n)
end

#=
Get intersection cuts
input:
- data::SBRPData is the instance
output:
- (cuts1, cuts2)::Tuple{Tuple{Arcs, Arcs}} is the tuple of arcs
=#
function getIntersectionCuts(data::SBRPData)::Tuple{Vector{Arcs},Vector{Arcs}}

    @debug "Intersection cuts separator"

    # setup
    @debug "Retrieving input"

    B::VVi = data.B
    A::Arcs = data.D.A
    cuts1::Vector{Arcs} = Vector{Arcs}()
    cuts2::Vector{Arcs} = Vector{Arcs}()
    cliques::Vector{VVi} = Vector{VVi}()

    # max clique model
    @debug "Creating MILP model"

    max_clique::Model = direct_model(CPLEX.Optimizer())
    #set_silent(max_clique)
    @variable(max_clique, z[block::Vi in B], Bin)
    @objective(max_clique, Max, sum(map(block::Vi -> z[block], B)))
    @constraint(max_clique,
        [(block, block′)::Tuple{Vi,Vi} in [(B[i], B[j]) for i::Int in 1:length(B) for j::Int in i+1:length(B) if isempty(∩(B[i], B[j]))]],
        1 >= z[block] + z[block′])
    @constraint(max_clique, sum(map(block::Vi -> z[block], B)) >= 2)

    # 0 = Balanced (Default)
    # 1 = Feasibility (Find integer solutions quickly)
    # 2 = Optimality (Focus on moving the Best Bound)
    # 3 = Best Bound (Aggressive theoretical bound improvement)
    # 4 = Hidden Feasibility (Aggressive search for difficult integer points)

    set_optimizer_attribute(max_clique, "CPXPARAM_Emphasis_MIP", 1)

    #
    iteration::Int = 1

    while true
        @debug "Running model"

        # set threads to 1
        _set_cplex_threads!(max_clique, 1)

        # solve
        optimize!(max_clique)

        # base case
        if termination_status(max_clique) == MOI.INFEASIBLE
            @debug "Infeasible model"
            break
        end

        # get clique
        @debug "Retain blocks"
        B′::VVi = filter(block::Vi -> value(z[block]) > 0.5, B)

        # check if it is a clique
        if all([!isempty(∩(block, block′)) for block′::Vi in B for block::Vi in B′ if block′ != block])
            error("It is not a clique")
        end

        # After the first optimal clique is found,
        # keep the solver from over-thinking the next one
        if iteration == 1
            set_optimizer_attribute(max_clique, "CPXPARAM_Preprocessing_Presolve", 0)
            # Use a very small MIP gap to stop early if finding the *absolute*
            # max clique becomes hard in later iterations
            #set_optimizer_attribute(max_clique, "CPXPARAM_MIP_Tolerances_MIPGap", 0.05)
        end

        # store clique
        push!(cliques, B′)

        # update clique model
        @constraint(max_clique, sum(map(block::Vi -> z[block], B′)) <= length(B′) - 1)

        # update iteration
        iteration += 1
    end

    # get nodes of each block
    Vb::Si = getBlocksNodes(data)
    nodes_blocks::Dict{Int,VVi} = Dict{Int,VVi}(map(i::Int -> i => filter(block::Vi -> i in block, B), collect(Vb)))

    # get cuts
    @debug "Iterate through cliques"

    for clique::VVi in cliques

        # get intersection
        @debug "Get clique intersection"

        intersection::Vi = ∩(clique...)

        ######## intersection cuts 1 ########
        @debug "Get intersection  cuts 1"

        isolated_intersection::Vi = setdiff(intersection, ∪(i for block::Vi in setdiff(B, clique) for i::Int in block))
        #    isolated_intersection = sediff(intersection, ∪(block... for block ∈ setdiff(B, clique)))
        if length(isolated_intersection) > 1
            push!(cuts1, Arcs([(i, j) for (i::Int, j::Int) in χ(isolated_intersection) if i != j]))
        end

        ######## intersection cuts 2 ########
        @debug "Get intersection  cuts 2"

        for block::Vi in clique
            # get arcs incident to nodes that belong exclusively to the block
            #      independent_arcs = [a for (i, blocks) in nodes_blocks for a in δ⁺(A, i) if ∧(length(blocks) == 1, block ∈ blocks)]
            covered_arcs = [a for i::Int in block for a::Arc in δ⁺(A, i) if ⊆(nodes_blocks[i], clique)]

            # edge case
            #      isempty(independent_arcs) && continue
            isempty(covered_arcs) && continue

            # store
            #      push!(cuts2, Arcs(∪([a for i ∈ intersection for a ∈ δ⁺(A, i)], independent_arcs)))
            push!(cuts2, Arcs(covered_arcs))
        end

    end

    return cuts1, cuts2
end

#=
Add intersection cuts type 1
input:
- model::Model is the MILP model
- cuts::Vector{Arcs} is the list of arcs
=#
function addIntersectionCuts1(model::Model, cuts::Vector{Arcs})
    for arcs::Arcs in cuts
        @constraint(model, sum(map(a::Arc -> model[:x][a], arcs)) == 0)
    end
end

#=
Add intersection cuts type 2
input:
- model::Model is the MILP model
- cuts::Vector{Arcs} is the list of arcs
=#
function addIntersectionCuts2(model::Model, cuts::Vector{Arcs})
    for arcs::Arcs in cuts
        @constraint(model, sum(map(a::Arc -> model[:x][a], arcs)) <= 2)
    end
end

#=
Add subtour inequalities
input:
- model::Model is a Mathematical Programming model
- sets::Set{Pair{Arcs, Arcs}} is a list of list of arcs
=#
function addSubtourCuts(model::Model, sets::Set{Tuple{Arcs,Arcs}})

    # get vars
    x = model[:x]

    @debug "Adding subtour cuts"

    for (lhs_arcs::Arcs, rhs_arcs::Arcs) in sets
        # add constraint
        @constraint(model, sum(lhs_arc::Arc -> x[lhs_arc], lhs_arcs) >= sum(rhs_arc::Arc -> x[rhs_arc], rhs_arcs) - 1)
    end
end

function addSubtourCut(model::Model, lhs_arcs::Arcs, rhs_arcs::Arcs)

    # get vars
    x = model[:x]

    # add constraint
    @constraint(model, sum(lhs_arc::Arc -> x[lhs_arc], lhs_arcs) >= sum(rhs_arc::Arc -> x[rhs_arc], rhs_arcs) - 1)

end

#=
Get subtour cuts
input:
- data::SBRPData is the instance
- model::Model is a Mathematical Programming model
- info::Dict{String, String} is the output log relation
output:
- components::Set{Tuple{Arcs, Arcs}} is the set of components to separate.
=#
function getSubtourCuts(
    data::SBRPData,
    model::Model,
    app::Dict{String,Any},
    info::Dict{String,String}
)::Set{Tuple{Arcs,Arcs}}

    @debug "Obtaining subtour cuts"

    x, z, w = model[:x], model[:z], model[:w]

    # setup
    V::Vi = collect(keys(data.D.V))
    A::Arcs = data.D.A
    B::VVi = data.B

    # outputs
    components::Set{Tuple{Arcs,Arcs}} = Set{Tuple{Arcs,Arcs}}()

    # helpers
    Vₘ = Dict{Int,Int}(map((idx, i)::Tuple{Int,Int} -> i => idx, enumerate(V)))
    Vₘʳ = Dict{Int,Int}(map((idx, i)::Tuple{Int,Int} -> idx => i, enumerate(V)))

    n::Int = length(Vₘ)
    iteration::Int = 1

    # Before the while loop
    set_optimizer_attribute(model, "CPXPARAM_LPMethod", 2)

    #
    set_optimizer_attribute(model, "CPXPARAM_Advance", 1)

    #
    epsilon::Float64 = 1e-2

    # get cuts greedly
    while true

        # Set silent mode
        #set_silent(model)

        # set threads to 1
        _set_cplex_threads!(model, 1)

        # store time
        time::Float64 = @elapsed optimize!(model)

        # error checking
        if !in(termination_status(model), [MOI.OPTIMAL, MOI.TIME_LIMIT, MOI.ALMOST_INFEASIBLE])
            throw(InvalidStateException("The model could not be solved"))
        end

        # get values
        w_val::Dict{Int,Float64} = Dict{Int,Float64}(map(i::Int -> i => value(w[i]), V))
        z_val::Dict{Int,Float64} = Dict{Int,Float64}(map(i::Int -> i => value(z[i]), V))
        x_val::ArcCostMap = ArcCostMap(map(a::Arc -> a => value(x[a]), A))

        # get subsets
        g = SparseMaxFlowMinCut.ArcFlow[]
        M::Int = 100000
        new_components::Set{Tuple{Arcs,Arcs}} = Set{Tuple{Arcs,Arcs}}()

        # used arcs
        A′::Arcs = filter(a::Arc -> x_val[a] > EPS, A)

        # used nodes
        V′::Vi = filter(i::Int -> z_val[i] > EPS, V)

        # depots
        depots′::Vi = filter(i::Int -> w_val[i] > EPS, V)

        # mounting graph
        for (i::Int, j::Int) in A′
            push!(g, SparseMaxFlowMinCut.ArcFlow(Vₘ[i], Vₘ[j], trunc(floor(x_val[Arc(i, j)], digits=5) * M)))
        end

        #
        max_violation::Float64 = 0.0

        # get subsets
        for source::Int in depots′

            for target::Int in V′

                # edge case
                source == target && continue

                # init
                maxFlow::Float64, flows, set = SparseMaxFlowMinCut.find_maxflow_mincut(
                    SparseMaxFlowMinCut.Graph(n, g),
                    Vₘ[source],
                    Vₘ[target]
                )
                flow::Float64 = maxFlow / M

                # base case: In the same component
                (set[Vₘ[target]] == 1) && continue

                # base case: Condition not met
                flow + epsilon >= z_val[source] + z_val[target] - 1 && continue

                # get set
                S::Si = Si(map(i::Int -> Vₘʳ[i], filter(i::Int -> set[i] == 1, 1:n)))

                # base case 3
                #      length(S) <= 1 && continue

                # get components
                Aₛ::Arcs = δ⁺(A, S)
                Aᵢ::Arcs = union(δ⁺(A′, source), δ⁺(A′, target))

                # calculate violation degree
                violation::Float64 = z_val[source] + z_val[target] - 1 - (flow + epsilon)


                # best improvement case
                if app["subcycle-separation"] == "best"

                    if max_violation < violation

                        # clean
                        empty!(new_components)

                    else

                        continue

                    end
                end

                # update
                max_violation = max(violation, max_violation)

                # store
                push!(new_components, (Aₛ, Aᵢ))
                addSubtourCut(model, Aₛ, Aᵢ)

                #
                #println("S: ", S)

                # first improvement case
                app["subcycle-separation"] == "first" && break

            end

            # first improvement case
            (app["subcycle-separation"] == "first" && !isempty(new_components)) && break

        end

        # base case
        isempty(new_components) && break

        # store components
        union!(components, new_components)

        # add ineqs
        #addSubtourCuts(model, components)

        # After the first iteration, you could potentially disable presolve
        # if the number of cuts is small and iterations are many.
        if iteration == 1
            set_optimizer_attribute(model, "CPXPARAM_Preprocessing_Presolve", 0)
            set_optimizer_attribute(model, "CPXPARAM_Preprocessing_Aggregator", 0) # Prevents re-combining rows
        end
        iteration += 1
    end

    #
    return components
end

"""Deep-copy intersection cut lists for caching across budget sweeps."""
function _copy_intersection_cuts(c1::Vector{Arcs}, c2::Vector{Arcs})::Tuple{Vector{Arcs},Vector{Arcs}}
    return (map(copy, c1), map(copy, c2))
end

"""
Set JuMP MIP start from a feasible `SBRPSolution` (hint only; may violate MTZ after `T` changes).
"""
function completeDigraphMipStart!(
    model::Model,
    data::SBRPData,
    sol::SBRPSolution,
    depot::Int,
    A::Arcs,
    V::Vi,
    B::VVi,
)::Nothing

    x = model[:x]
    y = model[:y]
    z = model[:z]
    w = model[:w]
    t = model[:t]
    Tlim::Float64 = data.T

    tour_nodes::Si = Set{Int}(sol.tour)
    block_node::Si = Set{Int}()
    for block::Vi in sol.B
        union!(block_node, block)
    end

    y_on(block::Vi)::Bool = any(sb::Vi -> sb == block, sol.B)

    for block::Vi in B
        set_start_value(y[block], y_on(block) ? 1.0 : 0.0)
    end

    for i::Int in V
        vis::Bool = (i in tour_nodes) || (i in block_node)
        set_start_value(z[i], vis ? 1.0 : 0.0)
    end

    start_depot::Int = isempty(sol.tour) ? depot : first(sol.tour)
    for i::Int in V
        set_start_value(w[i], i == start_depot ? 1.0 : 0.0)
    end

    tour_arc::Set{Arc} = Set{Arc}()
    for k::Int in 1:(length(sol.tour)-1)
        push!(tour_arc, sol.tour[k] => sol.tour[k+1])
    end

    for a::Arc in A
        set_start_value(x[a], (a in tour_arc) ? 1.0 : 0.0)
    end

    for i::Int in V
        set_start_value(t[i], 0.0)
    end
    set_start_value(t[depot], 0.0)

    cum::Float64 = 0.0
    for k::Int in 2:length(sol.tour)
        prev::Int = sol.tour[k-1]
        curr::Int = sol.tour[k]
        cum += time(data, prev => curr)
        set_start_value(t[curr], min(cum, Tlim))
    end

    return nothing
end

#=
Build SBRP model
input:
- data::SBRPData is a SBRP instance
- app::Dict{String, Any} is the relation of application parameters
=#
function runCOPCompleteDigraphIPModel(
    data::SBRPData,
    app::Dict{String,Any}
)::Tuple{SBRPSolution,Dict{String,String}}

    # instance parameters
    depot::Int = data.depot
    B::VVi = data.B
    A::Arcs = data.D.A
    A_set::ArcsSet = ArcsSet(data.D.A)
    T::Float64 = data.T
    V::Vi = collect(keys(data.D.V))
    profits::Dict{Vi,Float64} = data.profits
    Vb::Si = getBlocksNodes(data)
    info::Dict{String,String} = Dict{String,String}(
        "lazyCuts" => "0",
        "warmStartUsed" => "false",
        "pooledSubtourCutsSeeded" => "0",
        "usedCachedIntersection" => "n/a",
    )

    new_subtour_cuts::Set{Tuple{Arcs,Arcs}} = Set{Tuple{Arcs,Arcs}}()

    # new digraph
    A′::Arcs = filter((i, j)::Arc -> j != depot, A)
    V′::Si = setdiff(Si(V), depot)

    P′::Arcs = filter((i, j)::Arc -> i < j && Arc(j, i) in A_set, A)

    nodes_blocks::Dict{Int,VVi} = Dict{Int,VVi}(map(i::Int -> i => filter(block::Vi -> i in block, B), collect(Vb)))

    # Diagnostic: JuMP sum(x[a], δ±(A,i)) fails if either star is empty.
    empty_δm::Vector{Int} = Int[]
    empty_δp::Vector{Int} = Int[]
    for i::Int in V
        isempty(δ⁻(A, i)) && push!(empty_δm, i)
        isempty(δ⁺(A, i)) && push!(empty_δp, i)
    end
    if !isempty(empty_δm) || !isempty(empty_δp)
        n_show::Int = min(20, max(length(empty_δm), length(empty_δp)))
        @warn "runCOPCompleteDigraphIPModel: empty arc stars (will break degree/update_z constraints)" depot = depot n_V =
            length(V) n_A = length(A) n_empty_din = length(empty_δm) n_empty_dout = length(empty_δp) sample_din =
            empty_δm[1:min(end, n_show)] sample_dout = empty_δp[1:min(end, n_show)]
    else
        @debug "runCOPCompleteDigraphIPModel: every vertex has nonempty δ⁻ and δ⁺" depot = depot n_V = length(V) n_A = length(A)
    end

    model::Model = direct_model(CPLEX.Optimizer())
    unset_silent(model)

    cplex_time_limit_sec::Float64 = _ip_base_cplex_seconds(parse(Float64, get(app, "time-limit", "0")))
    set_time_limit_sec(model, cplex_time_limit_sec)

    @variable(model, x[a::Arc in A], Bin)
    @variable(model, y[b::Vi in B], Bin)
    @variable(model, z[i::Int in V], Bin)
    @variable(model, w[i::Int in V], Bin)
    @variable(model, t[i::Int in V], lower_bound = 0, upper_bound = T)

    @objective(model, Max, sum(block::Vi -> data.profits[block] * y[block], B))

    @constraint(model, sum(i::Int -> w[i], V) == 1)
    @constraint(model, lift_w[i::Int in V], w[i] <= z[i])
    @constraint(model, degree[i::Int in V], sum(a::Arc -> x[a], δ⁻(A, i)) == sum(a::Arc -> x[a], δ⁺(A, i)))
    @constraint(model, update_z[i::Int in V], sum(a::Arc -> x[a], δ⁻(A, i)) == z[i])
    @constraint(model, serviced_block[block::Vi in B], sum(a::Arc -> x[a], δ⁺(A, block)) >= y[block])
    @constraint(model, sum(a::Arc -> time(data, a) * x[a], A) <= T - sum(block::Vi -> y[block] * blockTime(data, block), B))

    # Nodes MTZ
    @constraint(model, t[depot] == 0.0)
    @constraint(model, mtz[a::Arc in A′], t[last(a)] >= t[first(a)] + x[a] * time(data, a) - (1 - x[a]) * T - x[reverse(a)] * time(data, reverse(a)))
    @constraint(model, ub1[i::Int in V′], t[i] <= T - sum(block::Vi -> y[block] * blockTime(data, block), B))
    @constraint(model, ub2[i::Int in V′], t[i] <= z[i] * T)

    # improvements
    @constraint(model, block3[block::Vi in B], y[block] - sum(x[a] for i::Int in block for a::Arc in δ⁺(A, i) if length(nodes_blocks[i]) == 1) >= 0)
    @constraint(model, subcycle_size_two[a::Arc in P′], x[a] + x[reverse(a)] <= 1)


    # getting intersection cuts
    if app["intersection-cuts"]
        @debug "Getting intersection cuts"

        ic_ref = get(app, "intersection_cuts_cache_ref", nothing)
        if haskey(app, "precomputed_intersection_cuts1") && haskey(app, "precomputed_intersection_cuts2")
            intersection_cuts1::Vector{Arcs} = app["precomputed_intersection_cuts1"]::Vector{Arcs}
            intersection_cuts2::Vector{Arcs} = app["precomputed_intersection_cuts2"]::Vector{Arcs}
            info["intersectionCutsTime"] = "cached"
            info["usedCachedIntersection"] = "true"
        else
            elapsed_time::Float64 = @elapsed begin
                intersection_cuts1, intersection_cuts2 = getIntersectionCuts(data)
            end
            info["intersectionCutsTime"] = string(elapsed_time)
            info["usedCachedIntersection"] = "false"
            if ic_ref !== nothing && ic_ref[] === nothing
                ic_ref[] = _copy_intersection_cuts(intersection_cuts1, intersection_cuts2)
            end
        end

        info["intersectionCuts1"], info["intersectionCuts2"] =
            string(length(intersection_cuts1)), string(length(intersection_cuts2))

        addIntersectionCuts1(model, intersection_cuts1)
        addIntersectionCuts2(model, intersection_cuts2)
    end

    # getting initial relaxation with both variables <x and y> relaxed
    @debug "Getting initial relaxation"

    # creating model with all variables relaxed
    unsetBinary(values(z))
    unsetBinary(values(w))
    unsetBinary(values(y))
    unsetBinary(values(x))

    _set_cplex_threads!(model, 1)

    optimize!(model)

    info["initialLP"] = string(objective_value(model))

    # getting initial relaxation with only x relaxed (y, w, and z integer)
    @debug "Getting initial relaxation with y, w, and z as integer"

    setBinary(values(z))
    setBinary(values(w))
    setBinary(values(y))

    info["yLPTime"] = string(@elapsed optimize!(model))
    info["yLP"] = string(objective_value(model))

    # get max-flow cuts with x and y relaxed or integer
    if app["subcycle-separation"] != "none"
        if !app["y-integer"]
            unsetBinary(values(y))
        end
        if !app["z-integer"]
            unsetBinary(values(z))
        end
        if !app["w-integer"]
            unsetBinary(values(w))
        end

        reuse_cuts::Bool = get(app, "reuse-cuts", false)
        pool::Union{Nothing,Set{Tuple{Arcs,Arcs}}} = get(app, "subtour_cut_pool", nothing)
        seed_pool::Bool = reuse_cuts && pool !== nothing && !isempty(pool)

        if seed_pool
            addSubtourCuts(model, pool::Set{Tuple{Arcs,Arcs}})
            info["pooledSubtourCutsSeeded"] = string(length(pool::Set{Tuple{Arcs,Arcs}}))
        else
            info["pooledSubtourCutsSeeded"] = "0"
        end

        info["maxFlowCutsTime"] = string(@elapsed new_subtour_cuts = getSubtourCuts(data, model, app, info))

        info["maxFlowCuts"] = string(length(new_subtour_cuts))

        # subtour cuts
        addSubtourCuts(model, new_subtour_cuts)
        optimize!(model)

        info["maxFlowLP"] = string(objective_value(model))
    else
        info["pooledSubtourCutsSeeded"] = "0"
    end

    # integer model
    setBinary(values(z))
    setBinary(values(w))
    setBinary(values(y))
    setBinary(values(x))

    warm_sol = get(app, "warm_start_solution", nothing)
    if warm_sol !== nothing
        completeDigraphMipStart!(model, data, warm_sol::SBRPSolution, depot, A, V, B)
        info["warmStartUsed"] = "true"
    else
        info["warmStartUsed"] = "false"
    end

    # run
    info["solverTime"] = string(@elapsed optimize!(model))

    # TIME_LIMIT (or other stop) with no integer incumbent → MOI has 0 results; objective_value throws.
    if !has_values(model)
        info["cost"] = "0.00"
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
        info["noFeasibleSolution"] = "true"
        pool_mut = get(app, "subtour_cut_pool", nothing)
        if get(app, "reuse-cuts", false) && pool_mut !== nothing
            union!(pool_mut::Set{Tuple{Arcs,Arcs}}, new_subtour_cuts)
        end
        return SBRPSolution(Vi([depot]), VVi()), info
    end

    info["cost"] = @sprintf("%.2f", objective_value(model))
    info["relativeGAP"] = string(relative_gap(model))
    info["nodeCount"] = string(node_count(model))

    # retrieve solution

    solution_arcs::Arcs = Arcs(filter(a::Arc -> value(x[a]) > 0.5, A))
    solution_nodes::Vi = Vi(filter(i::Int -> value(z[i]) > 0.5, V))
    solution_blocks::VVi = VVi(filter(block::Vi -> value(y[block]) > 0.5, B))
    chosen_depot::Int = depot

    tour::Vi = Vi([chosen_depot])

    while !isempty(solution_arcs)

        arc_idx::Union{Int,Nothing} = findfirst(a::Arc -> first(a) == last(tour), solution_arcs)

        if arc_idx == nothing
            error("Desired arc not found")
        end

        a::Arc = solution_arcs[arc_idx]

        i::Int, j::Int = first(a), last(a)

        push!(tour, j)

        deleteat!(solution_arcs, arc_idx)
    end

    solution::SBRPSolution = SBRPSolution(tour, solution_blocks)

    pool_mut = get(app, "subtour_cut_pool", nothing)
    if get(app, "reuse-cuts", false) && pool_mut !== nothing
        union!(pool_mut::Set{Tuple{Arcs,Arcs}}, new_subtour_cuts)
    end

    # return
    return solution, info
end
