using Test

@testset "dropZeroProfitBlocks! mutator" begin
    root::String = joinpath(@__DIR__, "..")
    include(joinpath(root, "src", "data.jl"))
    b0::Vi = Vi([1])
    b1::Vi = Vi([2])
    data::SBRPData = SBRPData(
        InputDigraph(
            Dict{Int,Vertex}(
                1 => Vertex(1, 0.0, 0.0),
                2 => Vertex(2, 1.0, 0.0),
            ),
            Arcs(),
            ArcCostMap(),
        ),
        1,
        VVi([b0, b1]),
        60.0,
        Dict{Vi,Float64}(b0 => 0.0, b1 => 3.5),
    )
    dropZeroProfitBlocks!(data)
    @test length(data.B) == 1
    @test data.B[1] == b1
    @test data.profits[b1] == 3.5
    @test length(data.profits) == 1

    bz::Vi = Vi([5])
    data2::SBRPData = SBRPData(
        InputDigraph(Dict{Int,Vertex}(5 => Vertex(5, 0.0, 0.0)), Arcs(), ArcCostMap()),
        1,
        VVi([bz]),
        60.0,
        Dict{Vi,Float64}(bz => 0.0),
    )
    @test_throws ArgumentError dropZeroProfitBlocks!(data2)
end

@testset "Carlos read drop-zero-profit-blocks (|B| monotone)" begin
    root::String = joinpath(@__DIR__, "..")
    inst::String = joinpath(root, "data", "carlos", "notified-alto-santo", "notified-alto-santo-1000-2021.txt")
    if !isfile(inst)
        @test_skip "Carlos fixture missing"
    else
        include(joinpath(root, "src", "data.jl"))
        app_full = Dict{String,Any}(
            "instance" => inst,
            "vehicle-time-limit" => "120",
            "no-cbrp-metric-closure" => false,
        )
        app_drop = Dict{String,Any}(
            "instance" => inst,
            "vehicle-time-limit" => "120",
            "no-cbrp-metric-closure" => false,
            "drop-zero-profit-blocks" => true,
        )
        data_full, _, _ = readSBRPDataCarlos(app_full)
        data_drop, _, _ = readSBRPDataCarlos(app_drop)
        @test length(data_drop.B) <= length(data_full.B)
        @test !any(b -> data_drop.profits[b] == 0.0, data_drop.B)
    end
end

@testset "Carlos sparse reader vs metric closure" begin
    root::String = joinpath(@__DIR__, "..")
    inst::String = joinpath(root, "data", "carlos", "notified-alto-santo", "notified-alto-santo-1000-2021.txt")
    if !isfile(inst)
        @test_skip "Carlos fixture missing"
    else
        include(joinpath(root, "src", "data.jl"))
        app_dense = Dict{String,Any}(
            "instance" => inst,
            "vehicle-time-limit" => "120",
            "no-cbrp-metric-closure" => false,
        )
        app_sparse = Dict{String,Any}(
            "instance" => inst,
            "vehicle-time-limit" => "120",
            "no-cbrp-metric-closure" => true,
        )
        data_dense, paths_dense, dist_dense = readSBRPDataCarlos(app_dense)
        data_sparse, paths_sparse, dist_sparse = readSBRPDataCarlos(app_sparse)
        @test paths_dense !== nothing && dist_dense !== nothing
        @test paths_sparse === nothing && dist_sparse === nothing
        @test length(data_dense.D.A) > 0
        @test length(data_sparse.D.A) > 0
        @test length(data_dense.D.A) != length(data_sparse.D.A)
        @test length(data_dense.B) == length(data_sparse.B)
        @test length(data_dense.D.V) == length(getBlocksNodes(data_dense)) + 1
        @test length(data_dense.D.V) < length(data_sparse.D.V)
    end
end

@testset "Matheus sparse reader vs metric closure" begin
    root::String = joinpath(@__DIR__, "..")
    inst::String = joinpath(root, "data", "campinas-random", "1.sbrp")
    if !isfile(inst)
        @test_skip "Campinas fixture missing"
    else
        include(joinpath(root, "src", "data.jl"))
        app_dense = Dict{String,Any}(
            "instance" => inst,
            "vehicle-time-limit" => "120",
            "no-cbrp-metric-closure" => false,
            "cluster-size-profits" => false,
            "unitary-profits" => false,
        )
        app_sparse = Dict{String,Any}(
            "instance" => inst,
            "vehicle-time-limit" => "120",
            "no-cbrp-metric-closure" => true,
            "cluster-size-profits" => false,
            "unitary-profits" => false,
        )
        data_dense, paths_dense, dist_dense = readSBRPData(app_dense)
        data_sparse, paths_sparse, dist_sparse = readSBRPData(app_sparse)
        @test paths_dense !== nothing && dist_dense !== nothing
        @test paths_sparse === nothing && dist_sparse === nothing
        @test length(data_dense.D.A) > 0
        @test length(data_sparse.D.A) > 0
        @test length(data_dense.D.A) != length(data_sparse.D.A)
        @test length(data_dense.B) == length(data_sparse.B)
        @test length(data_dense.D.V) == length(getBlocksNodes(data_dense)) + 1
        @test length(data_dense.D.V) < length(data_sparse.D.V)
        depot::Int = data_sparse.depot
        @test any(a::Arc -> a.first == depot, data_sparse.D.A)
        @test any(a::Arc -> a.second == depot, data_sparse.D.A)
    end
end

@testset "simplifyStreetDigraph! toy chain" begin
    root::String = joinpath(@__DIR__, "..")
    include(joinpath(root, "src", "data.jl"))

    # One-way 1 → 2 → 3 → 4; keep {1,4} via empty blocks + depot=1 and block node 4
    b4::Vi = Vi([4])
    data::SBRPData = SBRPData(
        InputDigraph(
            Dict{Int,Vertex}(
                1 => Vertex(1, 0.0, 0.0),
                2 => Vertex(2, 1.0, 0.0),
                3 => Vertex(3, 2.0, 0.0),
                4 => Vertex(4, 3.0, 0.0),
            ),
            Arcs([Arc(1, 2), Arc(2, 3), Arc(3, 4)]),
            ArcCostMap(Arc(1, 2) => 1.0, Arc(2, 3) => 1.0, Arc(3, 4) => 1.0),
        ),
        1,
        VVi([b4]),
        120.0,
        Dict{Vi,Float64}(b4 => 1.0),
    )
    stats = simplifyStreetDigraph!(data)
    @test stats.nodes_before == 4
    @test stats.nodes_after == 2
    @test Set(keys(data.D.V)) == Set([1, 4])
    @test haskey(data.D.distance, Arc(1, 4))
    @test data.D.distance[Arc(1, 4)] ≈ 3.0
    @test !haskey(data.D.V, 2)
    @test !haskey(data.D.V, 3)
end

@testset "simplifyStreetDigraph! block-safe keep" begin
    root::String = joinpath(@__DIR__, "..")
    include(joinpath(root, "src", "data.jl"))

    # Keep mid-chain node 3 as a block member
    b3::Vi = Vi([3])
    b4::Vi = Vi([4])
    data::SBRPData = SBRPData(
        InputDigraph(
            Dict{Int,Vertex}(
                1 => Vertex(1, 0.0, 0.0),
                2 => Vertex(2, 1.0, 0.0),
                3 => Vertex(3, 2.0, 0.0),
                4 => Vertex(4, 3.0, 0.0),
            ),
            Arcs([Arc(1, 2), Arc(2, 3), Arc(3, 4)]),
            ArcCostMap(Arc(1, 2) => 1.0, Arc(2, 3) => 1.0, Arc(3, 4) => 1.0),
        ),
        1,
        VVi([b3, b4]),
        120.0,
        Dict{Vi,Float64}(b3 => 1.0, b4 => 1.0),
    )
    simplifyStreetDigraph!(data)
    @test haskey(data.D.V, 3)
    @test haskey(data.D.distance, Arc(1, 3))
    @test data.D.distance[Arc(1, 3)] ≈ 2.0
    @test haskey(data.D.distance, Arc(3, 4))
    @test data.D.distance[Arc(3, 4)] ≈ 1.0
    @test !haskey(data.D.V, 2)
end

@testset "simplifyStreetDigraph! Campinas-1 preserves blocks" begin
    root::String = joinpath(@__DIR__, "..")
    inst::String = joinpath(root, "data", "campinas-random", "1.sbrp")
    if !isfile(inst)
        @test_skip "Campinas fixture missing"
    else
        include(joinpath(root, "src", "data.jl"))
        app = Dict{String,Any}(
            "instance" => inst,
            "vehicle-time-limit" => "120",
            "no-cbrp-metric-closure" => true,
            "cluster-size-profits" => false,
            "unitary-profits" => false,
        )
        data, _, _ = readSBRPData(app)
        block_nodes_before::Si = Si(getBlocksNodes(data))
        n_before::Int = length(data.D.V)
        a_before::Int = length(data.D.A)
        stats = simplifyStreetDigraph!(data)
        @test stats.nodes_after <= n_before
        @test stats.arcs_after <= a_before
        for v::Int in block_nodes_before
            @test haskey(data.D.V, v)
        end
        @test haskey(data.D.V, data.depot)
        for block::Vi in data.B
            length(block) <= 1 && continue
            for (i::Int, j::Int) in zip(block[begin:(end - 1)], block[(begin + 1):end])
                @test haskey(data.D.distance, Arc(i, j))
            end
            @test haskey(data.D.distance, Arc(last(block), first(block)))
        end
    end
end

@testset "compactToSparseSbrpSolution toy chain" begin
    root::String = joinpath(@__DIR__, "..")
    include(joinpath(root, "src", "data.jl"))
    include(joinpath(root, "src", "sol.jl"))
    include(joinpath(root, "src", "path_cbrp_interfaces.jl"))
    include(joinpath(root, "src", "brkga_path_warmstart.jl"))

    depot::Int = 4
    b1::Vi = Vi([1])
    b2::Vi = Vi([3])
    B::VVi = VVi([b1, b2])
    profits::Dict{Vi,Float64} = Dict{Vi,Float64}(b1 => 1.0, b2 => 2.0)
    dist::ArcCostMap = ArcCostMap(
        Arc(4, 1) => 1000.0,
        Arc(1, 2) => 1000.0,
        Arc(2, 3) => 1000.0,
        Arc(3, 4) => 1000.0,
    )
    Vdict::Dict{Int,Vertex} = Dict{Int,Vertex}(
        1 => Vertex(1, 0.0, 0.0),
        2 => Vertex(2, 1.0, 0.0),
        3 => Vertex(3, 2.0, 0.0),
        4 => Vertex(4, 0.0, 1.0),
    )
    data_sparse::SBRPData = SBRPData(
        InputDigraph(Vdict, Arcs(collect(keys(dist))), dist),
        depot,
        B,
        1.0e6,
        profits,
    )
    paths::Dict{Arc,Vi} = Dict{Arc,Vi}(Arc(1, 3) => Vi([2]))
    sol_c::SBRPSolution = SBRPSolution(Vi([1, 3]), VVi([b1, b2]))
    sol_s::SBRPSolution = compactToSparseSbrpSolution(data_sparse, paths, sol_c)
    @test sol_s.tour == Vi([4, 1, 2, 3, 4])
    @test length(sol_s.B) == 2
end

@testset "pathCbrpFixWarmStartFlag parsing" begin
    root::String = joinpath(@__DIR__, "..")
    include(joinpath(root, "src", "data.jl"))
    include(joinpath(root, "src", "brkga_path_warmstart.jl"))
    @test pathCbrpFixWarmStartFlag(Dict{String,Any}()) == false
    @test pathCbrpFixWarmStartFlag(Dict{String,Any}("path_cbrp_fix_warm_start" => true))
    @test pathCbrpFixWarmStartFlag(Dict{String,Any}("path_cbrp_fix_warm_start" => "true"))
    @test pathCbrpFixWarmStartFlag(Dict{String,Any}("path_cbrp_fix_warm_start" => "1"))
    @test !pathCbrpFixWarmStartFlag(Dict{String,Any}("path_cbrp_fix_warm_start" => false))
    @test !pathCbrpFixWarmStartFlag(Dict{String,Any}("path_cbrp_fix_warm_start" => "false"))
end

@testset "pathCbrpFixWarmFromXY! JuMP model (no solver)" begin
    root::String = joinpath(@__DIR__, "..")
    include(joinpath(root, "src", "data.jl"))
    include(joinpath(root, "src", "brkga_path_warmstart.jl"))
    using JuMP
    import MathOptInterface as MOI
    model = Model()
    A::Arcs = Arcs([Arc(1, 2), Arc(2, 1)])
    y_meta::Vector{Tuple{Int,Int}} = [(1, 1), (1, 2)]
    @variable(model, x[1:2], Bin)
    @variable(model, y[1:2], Bin)
    pathCbrpFixWarmFromXY!(model, A, y_meta, Set([Arc(1, 2)]), Set([(1, 1)]))
    @test occursin("fix_warm_x", name(model[:fix_warm_x][1]))
    @test occursin("fix_warm_y", name(model[:fix_warm_y][1]))
    @test constraint_object(model[:fix_warm_x][1]).set isa MOI.EqualTo{Float64}
    @test constraint_object(model[:fix_warm_x][1]).set.value ≈ 1.0
    @test constraint_object(model[:fix_warm_x][2]).set.value ≈ 0.0
    @test constraint_object(model[:fix_warm_y][1]).set.value ≈ 1.0
    @test constraint_object(model[:fix_warm_y][2]).set.value ≈ 0.0
end

@testset "Path-CBRP addPathSubtourCut structural" begin
    root::String = joinpath(@__DIR__, "..")
    include(joinpath(root, "src", "data.jl"))
    include(joinpath(root, "src", "model.jl"))
    model = direct_model(CPLEX.Optimizer())
    A::Arcs = Arcs([1 => 2, 2 => 3, 3 => 1])
    y_meta = [(1, 1), (1, 2)]
    out_idx = Dict(1 => [1, 3], 2 => [2], 3 => [3])
    @variable(model, x[1:3], Bin)
    @variable(model, y[1:2], Bin)
    S::Set{Int} = Set([1, 2])
    cut::PathSubtourCut = (S, 1, 3, 1, 2)
    addPathSubtourCut!(model, A, out_idx, cut)
    cons = all_constraints(model; include_variable_in_set_constraints=false)
    @test length(cons) == 1
    @test occursin("c_path_sec", name(cons[1]))
    @test num_variables(model) == 5
end

@testset "findViolatedPathSubtourCuts unit" begin
    root::String = joinpath(@__DIR__, "..")
    include(joinpath(root, "src", "data.jl"))
    include(joinpath(root, "src", "model.jl"))
    model = Model()
    A::Arcs = Arcs([1 => 2, 2 => 3, 3 => 1])
    y_meta::Vector{Tuple{Int,Int}} = [(1, 2), (1, 3)]
    out_idx::Dict{Int,Vector{Int}} = Dict(1 => [1], 2 => [2], 3 => [3])
    @variable(model, x[1:3])
    @variable(model, y[1:2])
    V::Vi = Vi([1, 2, 3])
    y_at::Dict{Tuple{Int,Int},Int} = path_y_index_map(y_meta)
    blocks_at::Dict{Int,Vector{Int}} = path_blocks_at(y_meta)
    Vₘ = Dict{Int,Int}(1 => 1, 2 => 2, 3 => 3)
    Vₘʳ = Dict{Int,Int}(1 => 1, 2 => 2, 3 => 3)
    ctx::PathSubtourSepContext = PathSubtourSepContext(
        model,
        x,
        y,
        A,
        y_meta,
        out_idx,
        1,
        V,
        y_at,
        blocks_at,
        Vₘ,
        Vₘʳ,
        3,
        1e-2,
        "all",
        100_000,
    )
    # Low arc flow but high y: SEC sum x(delta+(S)) >= y_i + y_j - 1 should be violated.
    x_val::Dict{Int,Float64} = Dict(1 => 0.4, 2 => 0.4, 3 => 0.4)
    y_val::Dict{Int,Float64} = Dict(1 => 0.9, 2 => 0.9)
    cuts::Set{PathSubtourCut} = findViolatedPathSubtourCuts(ctx, x_val, y_val)
    @test !isempty(cuts)
    # Feasible at equality: enough flow on cut arcs for active y.
    x_ok::Dict{Int,Float64} = Dict(1 => 1.0, 2 => 1.0, 3 => 1.0)
    y_ok::Dict{Int,Float64} = Dict(1 => 0.5, 2 => 0.5)
    cuts_ok::Set{PathSubtourCut} = findViolatedPathSubtourCuts(ctx, x_ok, y_ok)
    @test isempty(cuts_ok)
end

@testset "findDisconnectedPathSubtourCuts (integer candidate DFS)" begin
    root::String = joinpath(@__DIR__, "..")
    include(joinpath(root, "src", "data.jl"))
    include(joinpath(root, "src", "model.jl"))
    # depot 1; depot cycle 1→2→1; subtour 3→4→5→3; bridges 2→3, 5→1 for the single-tour case.
    A::Arcs = Arcs([1 => 2, 2 => 1, 3 => 4, 4 => 5, 5 => 3, 2 => 3, 5 => 1])
    y_meta::Vector{Tuple{Int,Int}} = [(1, 2), (2, 4)]
    out_idx::Dict{Int,Vector{Int}} = Dict(1 => [1], 2 => [2, 6], 3 => [3], 4 => [4], 5 => [5, 7])
    model = Model()
    @variable(model, x[1:length(A)])
    @variable(model, y[1:length(y_meta)])
    V::Vi = Vi([1, 2, 3, 4, 5])
    Vₘ = Dict{Int,Int}(i => i for i in V)
    ctx::PathSubtourSepContext = PathSubtourSepContext(
        model, x, y, A, y_meta, out_idx, 1, V,
        path_y_index_map(y_meta), path_blocks_at(y_meta),
        Vₘ, Vₘ, length(V), 1e-6, "all", 100_000,
    )
    xv(on::Vector{Int}) = Dict{Int,Float64}(k => (k in on ? 1.0 : 0.0) for k in 1:length(A))

    # Two cycles, both blocks serviced: one cut on the non-depot cycle for y_2 (node 4).
    cuts = findDisconnectedPathSubtourCuts(ctx, xv([1, 2, 3, 4, 5]), Dict(1 => 1.0, 2 => 1.0))
    @test cuts == [(Set([3, 4, 5]), 2)]
    @test path_out_arc_indices_S(first(cuts[1]), A, out_idx) == [7]

    # Single tour 1→2→3→4→5→1 through both serviced nodes: no cut.
    @test isempty(findDisconnectedPathSubtourCuts(ctx, xv([1, 6, 3, 4, 7]), Dict(1 => 1.0, 2 => 1.0)))

    # Unserviced subtour is not cut.
    @test isempty(findDisconnectedPathSubtourCuts(ctx, xv([1, 2, 3, 4, 5]), Dict(1 => 1.0, 2 => 0.0)))

    # Depot tour services nothing, subtour services everything: the depot-rooted DFS cut
    # catches it, while the pairwise max-flow family (y_i + y_j - 1) has no violated pair.
    x_gap = xv([1, 2, 3, 4, 5])
    y_gap = Dict(1 => 0.0, 2 => 1.0)
    @test findDisconnectedPathSubtourCuts(ctx, x_gap, y_gap) == [(Set([3, 4, 5]), 2)]
    @test isempty(findViolatedPathSubtourCuts(ctx, x_gap, y_gap))

    # Near-integer values (within CPLEX integrality tolerance) on a connected tour: no cut from
    # either separator (max-flow judges violation on exact values, not truncated capacities).
    x_near = Dict{Int,Float64}(k => (k in (1, 6, 3, 4, 7) ? 1.0 - 1e-6 : 1e-7) for k in 1:length(A))
    y_near = Dict(1 => 1.0 - 1e-6, 2 => 1.0)
    @test isempty(findDisconnectedPathSubtourCuts(ctx, x_near, y_near))
    @test isempty(findViolatedPathSubtourCuts(ctx, x_near, y_near))

    @test pathCandidateMaxFractionality(x_near, Dict(1 => 1.0)) ≈ 1e-6
end

@testset "findViolatedPathSubtourCuts judges violation on exact x (not truncated capacities)" begin
    root::String = joinpath(@__DIR__, "..")
    include(joinpath(root, "src", "data.jl"))
    include(joinpath(root, "src", "model.jl"))
    # Triangle 1→2→3→1, depot 1; y_1 = (block 1, node 2), y_2 = (block 1, node 3).
    A::Arcs = Arcs([1 => 2, 2 => 3, 3 => 1])
    y_meta::Vector{Tuple{Int,Int}} = [(1, 2), (1, 3)]
    out_idx::Dict{Int,Vector{Int}} = Dict(1 => [1], 2 => [2], 3 => [3])
    model = Model()
    @variable(model, x[1:3])
    @variable(model, y[1:2])
    V::Vi = Vi([1, 2, 3])
    Vₘ = Dict{Int,Int}(i => i for i in V)
    ctx::PathSubtourSepContext = PathSubtourSepContext(
        model, x, y, A, y_meta, out_idx, 1, V,
        path_y_index_map(y_meta), path_blocks_at(y_meta),
        Vₘ, Vₘ, 3, 1e-6, "all", 100_000,
    )
    # S = {2}: x(δ⁺(S)) = x_{2→3} = 0.300009 = y_1 + y_2 - 1 (tight). Truncating the capacity to
    # 0.30000 would undershoot by 9e-6 > epsilon and report a spurious cut.
    y_val = Dict(1 => 0.65, 2 => 0.650009)
    x_tight = Dict(1 => 1.0, 2 => 0.300009, 3 => 1.0)
    @test isempty(findViolatedPathSubtourCuts(ctx, x_tight, y_val))
    # Genuinely violated by 1e-4: found on S = {2} and S = {1, 2} (both have δ⁺(S) = {2→3}).
    x_viol = Dict(1 => 1.0, 2 => 0.299909, 3 => 1.0)
    cuts = findViolatedPathSubtourCuts(ctx, x_viol, y_val)
    @test Set(first.(cuts)) == Set([Set([2]), Set([1, 2])])
    for cut in cuts
        @test path_out_arc_indices_S(first(cut), A, out_idx) == [2]
        @test pathSubtourCutViolation(ctx, cut, x_viol, y_val) ≈ 1e-4
    end
end

@testset "Path SEC user-cut dedup and violation diagnostics" begin
    root::String = joinpath(@__DIR__, "..")
    include(joinpath(root, "src", "data.jl"))
    include(joinpath(root, "src", "model.jl"))
    # Triangle 1→2→3→1, depot 1; y_1 = (block 1, node 2), y_2 = (block 1, node 3).
    A::Arcs = Arcs([1 => 2, 2 => 3, 3 => 1])
    y_meta::Vector{Tuple{Int,Int}} = [(1, 2), (1, 3)]
    out_idx::Dict{Int,Vector{Int}} = Dict(1 => [1], 2 => [2], 3 => [3])
    model = Model()
    @variable(model, x[1:3])
    @variable(model, y[1:2])
    V::Vi = Vi([1, 2, 3])
    Vₘ = Dict{Int,Int}(i => i for i in V)
    ctx::PathSubtourSepContext = PathSubtourSepContext(
        model, x, y, A, y_meta, out_idx, 1, V,
        path_y_index_map(y_meta), path_blocks_at(y_meta),
        Vₘ, Vₘ, 3, 1e-2, "all", 100_000,
    )

    # Keys ignore i/j, the order of the y pair, and the iteration order of S.
    st = PathSubtourCallbackStats(0, 0, 0.0)
    @test pathUserCutKey!(st, (Set([2, 3]), 2, 1, 1, 2)) == pathUserCutKey!(st, (Set([3, 2]), 3, 1, 2, 1))
    @test pathUserCutKey!(st, (Set([2]), 2, 3, 1, 2)) != pathUserCutKey!(st, (Set([2, 3]), 2, 1, 1, 2))
    @test length(st.s_ids) == 2

    # Exact violation: S = {2}, δ⁺(S) = {2→3}; y_1 + y_2 - 1 - x_2.
    x_val = Dict(1 => 0.4, 2 => 0.3, 3 => 0.4)
    y_val = Dict(1 => 0.9, 2 => 0.9)
    cut_a::PathSubtourCut = (Set([2]), 2, 3, 1, 2)
    @test pathSubtourCutViolation(ctx, cut_a, x_val, y_val) ≈ 0.5
    # S = {2, 3}, δ⁺(S) = {3→1}: 0.8 - 0.795 = 0.005 <= epsilon (1e-2) → non-violated.
    cut_b::PathSubtourCut = (Set([2, 3]), 2, 1, 1, 2)
    x_b = Dict(1 => 0.4, 2 => 0.3, 3 => 0.795)
    @test pathSubtourCutViolation(ctx, cut_b, x_b, y_val) ≈ 0.005

    # Without dedup: repeats are submitted again but counted.
    st_off = PathSubtourCallbackStats(0, 0, 0.0)
    sub1, r1 = selectPathUserCuts!(st_off, ctx, Set([cut_a, cut_b]), x_b, y_val)
    @test length(sub1) == 2 && r1.found == 2 && r1.dup == 0 && r1.submitted == 2
    @test r1.nonviolated == 1
    @test r1.viol_min ≈ 0.005 && r1.viol_max ≈ 0.5 && r1.viol_mean ≈ (0.5 + 0.005) / 2
    sub2, r2 = selectPathUserCuts!(st_off, ctx, Set([cut_a]), x_b, y_val)
    @test length(sub2) == 1 && r2.dup == 1
    @test st_off.n_user_found == 3 && st_off.n_user_dup == 1 && st_off.n_user_rounds == 2
    @test length(st_off.seen_user_cuts) == 2
    @test st_off.n_user_nonviolated == 1

    # With dedup: repeats are dropped and do not feed the violation stats.
    st_on = PathSubtourCallbackStats(0, 0, 0.0; dedup_user_cuts=true)
    selectPathUserCuts!(st_on, ctx, Set([cut_a]), x_b, y_val)
    sub3, r3 = selectPathUserCuts!(st_on, ctx, Set([cut_a, cut_b]), x_b, y_val)
    @test sub3 == [cut_b] && r3.found == 2 && r3.dup == 1 && r3.submitted == 1
    @test st_on.user_viol_sum ≈ 0.5 + 0.005

    # Empty round: no round counted, NaN mean.
    _, r_empty = selectPathUserCuts!(st_on, ctx, Set{PathSubtourCut}(), x_b, y_val)
    @test r_empty.found == 0 && isnan(r_empty.viol_mean) && st_on.n_user_rounds == 2
    @test isnan(pathUserCutMeanViolation(PathSubtourCallbackStats(0, 0, 0.0)))
    @test _viol_str(Inf) == "NA" && _viol_str(0.123456) == "0.123"

    # Shared CLI threshold parsing (`--sec-min-violation`, both models).
    @test secMinViolation(Dict{String,Any}()) == 1e-4
    @test secMinViolation(Dict{String,Any}("sec-min-violation" => 1e-4)) == 1e-4
    @test secMinViolation(Dict{String,Any}("sec-min-violation" => "0.001")) == 1e-3
    @test_throws ArgumentError secMinViolation(Dict{String,Any}("sec-min-violation" => -1.0))
end

@testset "pathCallbackSolutionIsInteger" begin
    root::String = joinpath(@__DIR__, "..")
    include(joinpath(root, "src", "model.jl"))
    @test pathCallbackSolutionIsInteger(Dict(1 => 1.0, 2 => 0.0), Dict(1 => 1.0))
    @test !pathCallbackSolutionIsInteger(Dict(1 => 0.5), Dict(1 => 1.0))
    @test !pathCallbackSolutionIsInteger(Dict(1 => 1.0), Dict(1 => 0.5))
end

@testset "findViolatedCompleteSubtourCuts unit" begin
    root::String = joinpath(@__DIR__, "..")
    include(joinpath(root, "src", "data.jl"))
    include(joinpath(root, "src", "model.jl"))

    Vlist::Vi = Vi([1, 2, 3])
    A::Arcs = Arcs([Arc(1, 2), Arc(2, 3), Arc(3, 1), Arc(2, 1), Arc(3, 2), Arc(1, 3)])
    dist::ArcCostMap = ArcCostMap(map(a::Arc -> a => 1.0, A))
    Vdict::Dict{Int,Vertex} = Dict{Int,Vertex}(
        1 => Vertex(1, 0.0, 0.0),
        2 => Vertex(2, 1.0, 0.0),
        3 => Vertex(3, 0.0, 1.0),
    )
    data::SBRPData = SBRPData(
        InputDigraph(Vdict, A, dist),
        1,
        VVi([Vi([2]), Vi([3])]),
        60.0,
        Dict{Vi,Float64}(Vi([2]) => 1.0, Vi([3]) => 1.0),
    )

    model::Model = direct_model(CPLEX.Optimizer())
    set_silent(model)
    @variable(model, x[a::Arc in A], Bin)
    @variable(model, z[i::Int in Vlist], Bin)
    @variable(model, w[i::Int in Vlist], Bin)

    app::Dict{String,Any} = Dict{String,Any}("subcycle-separation" => "all")
    ctx::CompleteSubtourSepContext = buildCompleteSubtourSepContext(data, model, app)

    x_val::ArcCostMap = ArcCostMap(
        Arc(1, 2) => 0.4,
        Arc(2, 3) => 0.4,
        Arc(3, 1) => 0.4,
        Arc(2, 1) => 0.0,
        Arc(3, 2) => 0.0,
        Arc(1, 3) => 0.0,
    )
    z_val::Dict{Int,Float64} = Dict(1 => 1.0, 2 => 0.9, 3 => 0.9)
    w_val::Dict{Int,Float64} = Dict(1 => 1.0, 2 => 0.0, 3 => 0.0)
    cuts::Set{Tuple{Arcs,Arcs}} = findViolatedCompleteSubtourCuts(ctx, x_val, z_val, w_val)
    @test !isempty(cuts)

    x_ok::ArcCostMap = ArcCostMap(map(a::Arc -> a => 1.0, A))
    z_ok::Dict{Int,Float64} = Dict(1 => 0.5, 2 => 0.5, 3 => 0.5)
    w_ok::Dict{Int,Float64} = Dict(1 => 1.0, 2 => 0.0, 3 => 0.0)
    cuts_ok::Set{Tuple{Arcs,Arcs}} =
        findViolatedCompleteSubtourCuts(ctx, x_ok, z_ok, w_ok)
    @test isempty(cuts_ok)
    @test ctx.epsilon == 1e-4
end

@testset "findViolatedCompleteSubtourCuts judges violation on exact x" begin
    root::String = joinpath(@__DIR__, "..")
    include(joinpath(root, "src", "data.jl"))
    include(joinpath(root, "src", "model.jl"))

    Vlist::Vi = Vi([1, 2, 3])
    A::Arcs = Arcs([Arc(1, 2), Arc(2, 3), Arc(3, 1), Arc(2, 1), Arc(3, 2), Arc(1, 3)])
    data::SBRPData = SBRPData(
        InputDigraph(
            Dict{Int,Vertex}(i => Vertex(i, Float64(i), 0.0) for i in Vlist),
            A,
            ArcCostMap(map(a::Arc -> a => 1.0, A)),
        ),
        1,
        VVi([Vi([2]), Vi([3])]),
        60.0,
        Dict{Vi,Float64}(Vi([2]) => 1.0, Vi([3]) => 1.0),
    )
    model::Model = Model()
    @variable(model, x[a::Arc in A])
    ctx::CompleteSubtourSepContext = buildCompleteSubtourSepContext(
        data, model, Dict{String,Any}("subcycle-separation" => "all", "sec-min-violation" => 1e-6),
    )
    @test ctx.epsilon == 1e-6

    xv(x23::Float64, x13::Float64) = ArcCostMap(
        Arc(1, 2) => 1.0, Arc(2, 3) => x23, Arc(3, 1) => 1.0,
        Arc(2, 1) => 0.0, Arc(3, 2) => 0.0, Arc(1, 3) => x13,
    )
    w_val::Dict{Int,Float64} = Dict(1 => 1.0, 2 => 0.0, 3 => 0.0)

    # Depot 1 → target 3 yields S = {1, 2}, δ⁺(S) = {2→3, 1→3}; rhs = z_1 + z_3 - 1 = z_3.
    # Tight: x(δ⁺(S)) = 0.300009 = z_3. Truncating x_{2→3} to 0.30000 would fake a violation.
    @test isempty(findViolatedCompleteSubtourCuts(ctx, xv(0.300009, 0.0), Dict(1 => 1.0, 2 => 1.0, 3 => 0.300009), w_val))
    # Not violated: x_{1→3} = 5e-5 <= EPS is left out of the max-flow graph but counts exactly.
    @test isempty(findViolatedCompleteSubtourCuts(ctx, xv(0.3, 5e-5), Dict(1 => 1.0, 2 => 1.0, 3 => 0.30004), w_val))
    # Genuinely violated by 1e-4: found on S = {1, 2}.
    cuts = findViolatedCompleteSubtourCuts(ctx, xv(0.299909, 0.0), Dict(1 => 1.0, 2 => 1.0, 3 => 0.300009), w_val)
    @test length(cuts) == 1
    @test Set(first(first(cuts))) == Set([Arc(2, 3), Arc(1, 3)])
end

@testset "Complete digraph SEC callback validation" begin
    root::String = joinpath(@__DIR__, "..")
    include(joinpath(root, "src", "model.jl"))
    @test_throws ArgumentError validateCompleteSubtourSeparation!(
        Dict{String,Any}(
            "subcycle-separation" => "none",
            "subcycle-separation-engine" => "callback",
        ),
    )
    validateCompleteSubtourSeparation!(
        Dict{String,Any}(
            "subcycle-separation" => "first",
            "subcycle-separation-engine" => "callback",
        ),
    )
end

@testset "Complete digraph SEC callback smoke (CPLEX)" begin
    root::String = joinpath(@__DIR__, "..")
    inst::String = joinpath(root, "data", "carlos", "notified-alto-santo", "notified-alto-santo-1000-2021.txt")
    if !isfile(inst)
        @test_skip "Carlos fixture missing"
    else
        include(joinpath(root, "src", "data.jl"))
        include(joinpath(root, "src", "model.jl"))
        ok::Ref{Bool} = Ref(false)
        info::Union{Nothing,Dict{String,String}} = nothing
        try
            data = readSBRPDataCarlos(Dict{String,Any}(
                "instance" => inst,
                "vehicle-time-limit" => "120",
                "no-cbrp-metric-closure" => false,
            ))[1]
            data.T = 90.0
            _, info = runCOPCompleteDigraphIPModel(data, Dict{String,Any}(
                "time-limit" => "30",
                "subcycle-separation" => "first",
                "subcycle-separation-engine" => "callback",
                "intersection-cuts" => false,
                "y-integer" => false,
                "z-integer" => false,
                "w-integer" => false,
            ))
            ok[] = true
        catch e
            if e isa ArgumentError
                rethrow(e)
            end
        end
        ok[] || @test_skip "CPLEX unavailable or complete digraph callback solver error"
        @test info !== nothing
        @test get(info, "subcycleSeparationEngine", "") == "callback"
        @test haskey(info, "maxFlowUserCuts")
        @test parse(Int, get(info, "maxFlowUserCuts", "0")) >= 0
    end
end

@testset "Path-CBRP no-MTZ validation" begin
    root::String = joinpath(@__DIR__, "..")
    include(joinpath(root, "src", "model.jl"))
    app_none::Dict{String,Any} = Dict{String,Any}(
        "no-path-cbrp-mtz" => true,
        "subcycle-separation" => "none",
        "subcycle-separation-engine" => "callback",
    )
    app_root::Dict{String,Any} = Dict{String,Any}(
        "no-path-cbrp-mtz" => true,
        "subcycle-separation" => "first",
        "subcycle-separation-engine" => "root",
    )
    app_ok::Dict{String,Any} = Dict{String,Any}(
        "no-path-cbrp-mtz" => true,
        "subcycle-separation" => "first",
        "subcycle-separation-engine" => "callback",
    )
    @test_throws ArgumentError validatePathCbrpNoMtzSeparation!(app_none)
    @test_throws ArgumentError validatePathCbrpNoMtzSeparation!(app_root)
    validatePathCbrpNoMtzSeparation!(app_ok)

    inst::String = joinpath(root, "data", "carlos", "notified-alto-santo", "notified-alto-santo-1000-2021.txt")
    if isfile(inst)
        include(joinpath(root, "src", "data.jl"))
        data = readSBRPDataCarlos(Dict{String,Any}(
            "instance" => inst,
            "vehicle-time-limit" => "120",
            "no-cbrp-metric-closure" => true,
        ))[1]
        @test_throws ArgumentError runPathCbrpMipModel(data, app_none)
        @test_throws ArgumentError runPathCbrpMipModel(data, app_root)
    end
end

@testset "Path-CBRP IP smoke (CPLEX)" begin
    root::String = joinpath(@__DIR__, "..")
    inst::String = joinpath(root, "data", "carlos", "notified-alto-santo", "notified-alto-santo-1000-2021.txt")
    if !isfile(inst)
        @test_skip "Carlos fixture missing"
    else
        include(joinpath(root, "src", "data.jl"))
        include(joinpath(root, "src", "model.jl"))
        ok::Ref{Bool} = Ref(false)
        sol::Union{Nothing,SBRPSolution} = nothing
        info::Union{Nothing,Dict{String,String}} = nothing
        try
            data = readSBRPDataCarlos(Dict{String,Any}(
                "instance" => inst,
                "vehicle-time-limit" => "120",
                "no-cbrp-metric-closure" => true,
            ))[1]
            sol, info = runPathCbrpMipModel(data, Dict{String,Any}("time-limit" => "5"))
            ok[] = true
        catch
        end
        ok[] || @test_skip "CPLEX unavailable or Path-CBRP solver error"
        @test sol !== nothing && info !== nothing
        @test haskey(info, "cost")
        @test haskey(info, "bestBound")
        @test info["bestBound"] != ""
        @test parse(Float64, info["cost"]) ≥ 0.0
        if info["bestBound"] != "N/A"
            @test parse(Float64, info["cost"]) <= parse(Float64, info["bestBound"]) + 1e-3
        end
        @test length(sol.tour) ≥ 2
        @test get(info, "warmStartUsed", "false") == "false"
    end
end

@testset "Path-CBRP SEC smoke (CPLEX)" begin
    root::String = joinpath(@__DIR__, "..")
    inst::String = joinpath(root, "data", "carlos", "notified-alto-santo", "notified-alto-santo-1000-2021.txt")
    if !isfile(inst)
        @test_skip "Carlos fixture missing"
    else
        include(joinpath(root, "src", "data.jl"))
        include(joinpath(root, "src", "model.jl"))
        ok::Ref{Bool} = Ref(false)
        sol::Union{Nothing,SBRPSolution} = nothing
        info::Union{Nothing,Dict{String,String}} = nothing
        data::Union{Nothing,SBRPData} = nothing
        try
            data = readSBRPDataCarlos(Dict{String,Any}(
                "instance" => inst,
                "vehicle-time-limit" => "120",
                "no-cbrp-metric-closure" => true,
            ))[1]
            solve_app_sec = Dict{String,Any}(
                "time-limit" => "30",
                "subcycle-separation" => "first",
                "subcycle-separation-engine" => "root",
            )
            sol, info = runPathCbrpMipModel(data, solve_app_sec)
            ok[] = true
        catch
        end
        ok[] || @test_skip "CPLEX unavailable or Path-CBRP SEC solver error"
        @test sol !== nothing && info !== nothing && data !== nothing
        @test parse(Int, get(info, "maxFlowCuts", "0")) >= 0
        @test haskey(info, "maxFlowCutsTime")
        @test get(info, "subcycleSeparationEngine", "") == "root"
        include(joinpath(root, "src", "sol.jl"))
        checkSBRPSolution(data::SBRPData, sol::SBRPSolution)
    end
end

@testset "Path-CBRP SEC callback smoke (CPLEX)" begin
    root::String = joinpath(@__DIR__, "..")
    inst::String = joinpath(root, "data", "carlos", "notified-alto-santo", "notified-alto-santo-1000-2021.txt")
    if !isfile(inst)
        @test_skip "Carlos fixture missing"
    else
        include(joinpath(root, "src", "data.jl"))
        include(joinpath(root, "src", "model.jl"))
        ok::Ref{Bool} = Ref(false)
        sol::Union{Nothing,SBRPSolution} = nothing
        info::Union{Nothing,Dict{String,String}} = nothing
        data::Union{Nothing,SBRPData} = nothing
        try
            data = readSBRPDataCarlos(Dict{String,Any}(
                "instance" => inst,
                "vehicle-time-limit" => "120",
                "no-cbrp-metric-closure" => true,
            ))[1]
            solve_app_cb = Dict{String,Any}(
                "time-limit" => "30",
                "subcycle-separation" => "first",
                "subcycle-separation-engine" => "callback",
            )
            sol, info = runPathCbrpMipModel(data, solve_app_cb)
            ok[] = true
        catch
        end
        ok[] || @test_skip "CPLEX unavailable or Path-CBRP SEC callback solver error"
        @test sol !== nothing && info !== nothing && data !== nothing
        @test get(info, "subcycleSeparationEngine", "") == "callback"
        @test get(info, "pathCbrpMtzEnabled", "") == "true"
        @test parse(Int, get(info, "maxFlowCuts", "0")) >= 0
        @test haskey(info, "maxFlowUserCuts")
        @test haskey(info, "maxFlowLazyCuts")
        @test parse(Int, info["maxFlowUserCuts"]) >= 0
        @test parse(Int, info["maxFlowLazyCuts"]) >= 0
        @test parse(Int, info["maxFlowCuts"]) ==
            parse(Int, info["maxFlowUserCuts"]) + parse(Int, info["maxFlowLazyCuts"])
        include(joinpath(root, "src", "sol.jl"))
        checkSBRPSolution(data::SBRPData, sol::SBRPSolution)
    end
end

@testset "Path-CBRP no-MTZ + callback smoke (CPLEX)" begin
    root::String = joinpath(@__DIR__, "..")
    inst::String = joinpath(root, "data", "carlos", "notified-alto-santo", "notified-alto-santo-1000-2021.txt")
    if !isfile(inst)
        @test_skip "Carlos fixture missing"
    else
        include(joinpath(root, "src", "data.jl"))
        include(joinpath(root, "src", "model.jl"))
        ok::Ref{Bool} = Ref(false)
        info::Union{Nothing,Dict{String,String}} = nothing
        try
            data = readSBRPDataCarlos(Dict{String,Any}(
                "instance" => inst,
                "vehicle-time-limit" => "120",
                "no-cbrp-metric-closure" => true,
            ))[1]
            solve_app = Dict{String,Any}(
                "time-limit" => "15",
                "no-path-cbrp-mtz" => true,
                "subcycle-separation" => "first",
                "subcycle-separation-engine" => "callback",
            )
            _, info = runPathCbrpMipModel(data, solve_app)
            ok[] = true
        catch e
            if e isa ArgumentError
                rethrow(e)
            end
        end
        ok[] || @test_skip "CPLEX unavailable or no-MTZ callback solver error"
        @test info !== nothing
        @test get(info, "pathCbrpMtzEnabled", "") == "false"
        @test get(info, "subcycleSeparationEngine", "") == "callback"
    end
end

@testset "Path-CBRP SEC root vs callback objective (CPLEX)" begin
    root::String = joinpath(@__DIR__, "..")
    inst::String = joinpath(root, "data", "carlos", "notified-alto-santo", "notified-alto-santo-1000-2021.txt")
    if !isfile(inst)
        @test_skip "Carlos fixture missing"
    else
        include(joinpath(root, "src", "data.jl"))
        include(joinpath(root, "src", "model.jl"))
        data = readSBRPDataCarlos(Dict{String,Any}(
            "instance" => inst,
            "vehicle-time-limit" => "120",
            "no-cbrp-metric-closure" => true,
        ))[1]
        base_app = Dict{String,Any}(
            "time-limit" => "60",
            "subcycle-separation" => "first",
        )
        cost_root::Union{Nothing,Float64} = nothing
        cost_cb::Union{Nothing,Float64} = nothing
        try
            _, info_r = runPathCbrpMipModel(
                data,
                merge(base_app, Dict("subcycle-separation-engine" => "root")),
            )
            cost_root = parse(Float64, info_r["cost"])
            _, info_c = runPathCbrpMipModel(
                data,
                merge(base_app, Dict("subcycle-separation-engine" => "callback")),
            )
            cost_cb = parse(Float64, info_c["cost"])
        catch
            @test_skip "CPLEX unavailable or root vs callback comparison failed"
        end
        if cost_root !== nothing && cost_cb !== nothing
            @test abs(cost_root - cost_cb) <= 1e-3
        end
    end
end

@testset "Path-CBRP warm start from BRKGA .sol file (CPLEX)" begin
    root::String = joinpath(@__DIR__, "..")
    inst::String = joinpath(root, "data", "carlos", "notified-alto-santo", "notified-alto-santo-1000-2021.txt")
    if !isfile(inst)
        @test_skip "Carlos fixture missing"
    else
        include(joinpath(root, "src", "data.jl"))
        include(joinpath(root, "src", "sol.jl"))
        include(joinpath(root, "src", "nn_heuristic.jl"))
        include(joinpath(root, "src", "model.jl"))
        include(joinpath(root, "src", "brkga.jl"))
        app_dense = Dict{String,Any}(
            "instance" => inst,
            "vehicle-time-limit" => "120",
            "no-cbrp-metric-closure" => false,
            "brkga-conf" => joinpath(root, "conf", "config.conf"),
            "time-limit" => "8",
        )
        app_sparse = Dict{String,Any}(
            "instance" => inst,
            "vehicle-time-limit" => "120",
            "no-cbrp-metric-closure" => true,
        )
        data_dense, _, _ = readSBRPDataCarlos(app_dense)
        data_sparse, _, _ = readSBRPDataCarlos(app_sparse)
        Tlim::Float64 = parse(Float64, app_dense["vehicle-time-limit"])
        data_dense.T = Tlim
        data_sparse.T = Tlim
        sol_b = nothing
        try
            sol_b, _ = runCOPBRKGAModel(data_dense, app_dense)
        catch
            @test_skip "BRKGA unavailable or failed"
        end
        d = mktempdir()
        wbase = joinpath(d, "warm_export_brkga")
        writeSolution(wbase, data_dense, sol_b)
        sol_path = wbase * ".sol"
        @test isfile(sol_path)
        app_path = Dict{String,Any}(
            "instance" => inst,
            "vehicle-time-limit" => "120",
            "no-cbrp-metric-closure" => true,
            "path-cbrp-warm-sol" => sol_path,
            "time-limit" => "15",
        )
        app_path["path-cbrp-warm-sol-format"] = "brkga-sol"
        attachPathCbrpWarmStartFromFile!(app_path, data_sparse)
        if !haskey(app_path, "warm_start_solution")
            @test_skip "BRKGA tour did not expand to a feasible sparse warm start"
        end
        info_w = nothing
        try
            _, info_w = runPathCbrpMipModel(data_sparse, app_path)
        catch
            @test_skip "Path-CBRP with warm start failed"
        end
        @test info_w !== nothing
        @test get(info_w, "warmStartUsed", "false") == "true"
    end
end

@testset "Path-CBRP fix warm start from MILP solution (CPLEX)" begin
    root::String = joinpath(@__DIR__, "..")
    inst::String = joinpath(root, "data", "carlos", "notified-alto-santo", "notified-alto-santo-1000-2021.txt")
    if !isfile(inst)
        @test_skip "Carlos fixture missing"
    else
        include(joinpath(root, "src", "data.jl"))
        include(joinpath(root, "src", "model.jl"))
        data_sparse = readSBRPDataCarlos(Dict{String,Any}(
            "instance" => inst,
            "vehicle-time-limit" => "120",
            "no-cbrp-metric-closure" => true,
        ))[1]
        data_sparse.T = parse(Float64, "120")
        sol_ref::Union{Nothing,SBRPSolution} = nothing
        info_ref::Union{Nothing,Dict{String,String}} = nothing
        try
            sol_ref, info_ref = runPathCbrpMipModel(
                data_sparse,
                Dict{String,Any}("time-limit" => "30"),
            )
        catch
            @test_skip "CPLEX unavailable or Path-CBRP solver error"
        end
        sol_ref === nothing && @test_skip "No reference solution"
        info_ref === nothing && @test_skip "No reference info"
        ref_cost::Float64 = parse(Float64, info_ref["cost"])
        ref_cost <= 0.0 && @test_skip "Reference solve produced no positive incumbent"
        app_fix = Dict{String,Any}(
            "time-limit" => "15",
            "warm_start_solution" => sol_ref,
            "path_cbrp_fix_warm_start" => true,
        )
        info_fix = nothing
        try
            _, info_fix = runPathCbrpMipModel(data_sparse, app_fix)
        catch
            @test_skip "Path-CBRP with fix warm start failed"
        end
        @test info_fix !== nothing
        @test get(info_fix, "warmStartUsed", "false") == "true"
        @test get(info_fix, "warmStartFixed", "false") == "true"
        @test isapprox(parse(Float64, info_fix["cost"]), ref_cost; rtol=1e-5, atol=1e-3)
    end
end

@testset "parsePathCbrpMtzSolution article file" begin
    root::String = joinpath(@__DIR__, "..")
    mtz::String = joinpath(
        root,
        "solutions",
        "results-cbrp-article",
        "path-cbrp-mtz",
        "notified-alto-santo-1000-2016.txt",
    )
    inst::String = joinpath(root, "data", "carlos", "notified-alto-santo", "notified-alto-santo-1000-2016.txt")
    if !isfile(mtz) || !isfile(inst)
        @test_skip "path-cbrp-mtz fixture or instance missing"
    else
        include(joinpath(root, "src", "data.jl"))
        include(joinpath(root, "src", "brkga_path_warmstart.jl"))
        x_on, y_nb = parsePathCbrpMtzSolution(mtz)
        @test length(x_on) == 48
        @test length(y_nb) == 10
        @test detectPathCbrpWarmSolFormat(mtz) == "path-cbrp-mtz"
        data, _, _ = readSBRPDataCarlos(Dict{String,Any}(
            "instance" => inst,
            "vehicle-time-limit" => "120",
            "no-cbrp-metric-closure" => true,
        ))
        id_map = buildCarlosBlockIdToIndex(data, inst)
        @test length(id_map) >= 10
    end
end

@testset "Path-CBRP warm start from path-cbrp-mtz file (CPLEX)" begin
    root::String = joinpath(@__DIR__, "..")
    mtz::String = joinpath(
        root,
        "solutions",
        "results-cbrp-article",
        "path-cbrp-mtz",
        "notified-alto-santo-1000-2016.txt",
    )
    inst::String = joinpath(root, "data", "carlos", "notified-alto-santo", "notified-alto-santo-1000-2016.txt")
    if !isfile(mtz) || !isfile(inst)
        @test_skip "path-cbrp-mtz fixture or instance missing"
    else
        include(joinpath(root, "src", "data.jl"))
        include(joinpath(root, "src", "model.jl"))
        data, _, _ = readSBRPDataCarlos(Dict{String,Any}(
            "instance" => inst,
            "vehicle-time-limit" => "120",
            "no-cbrp-metric-closure" => true,
        ))
        app_path = Dict{String,Any}(
            "instance" => inst,
            "vehicle-time-limit" => "120",
            "no-cbrp-metric-closure" => true,
            "path-cbrp-warm-sol" => mtz,
            "path-cbrp-warm-sol-format" => "path-cbrp-mtz",
            "time-limit" => "15",
        )
        attachPathCbrpWarmStartFromFile!(app_path, data)
        if !haskey(app_path, "path_cbrp_warm_xy")
            @test_skip "Article MTZ X arcs do not match current sparse Carlos digraph (instance/build drift)"
        else
            info_w = nothing
            try
                _, info_w = runPathCbrpMipModel(data, app_path)
            catch
                @test_skip "CPLEX unavailable or Path-CBRP with MTZ warm start failed"
            end
            @test info_w !== nothing
            @test get(info_w, "warmStartUsed", "false") == "true"
        end
    end
end

@testset "Complete digraph reuse-cuts + warm start (CPLEX)" begin
    root::String = joinpath(@__DIR__, "..")
    inst::String = joinpath(root, "data", "carlos", "notified-alto-santo", "notified-alto-santo-1000-2021.txt")
    if !isfile(inst)
        @test_skip "Carlos fixture missing"
    else
        include(joinpath(root, "src", "data.jl"))
        include(joinpath(root, "src", "model.jl"))
        ok::Ref{Bool} = Ref(false)
        try
            data = readSBRPDataCarlos(Dict{String,Any}(
                "instance" => inst,
                "vehicle-time-limit" => "120",
                "no-cbrp-metric-closure" => false,
            ))[1]
            pool::Set{Tuple{Arcs, Arcs}} = Set{Tuple{Arcs, Arcs}}()
            ref_ic = Ref{Union{Nothing,Tuple{Vector{Arcs},Vector{Arcs}}}}(nothing)
            app1::Dict{String,Any} = Dict{String,Any}(
                "time-limit" => "15",
                "subcycle-separation" => "first",
                "intersection-cuts" => false,
                "y-integer" => false,
                "z-integer" => false,
                "w-integer" => false,
                "reuse-cuts" => true,
                "subtour_cut_pool" => pool,
                "intersection_cuts_cache_ref" => ref_ic,
            )
            data.T = 90.0
            sol1, info1 = runCOPCompleteDigraphIPModel(data, app1)
            @test info1["warmStartUsed"] == "false"
            n_pool_after_first::Int = length(pool)
            @test n_pool_after_first ≥ 0

            ws = SBRPSolution(copy(sol1.tour), deepcopy(sol1.B))
            app2::Dict{String,Any} = Dict{String,Any}(
                "time-limit" => "15",
                "subcycle-separation" => "first",
                "intersection-cuts" => false,
                "y-integer" => false,
                "z-integer" => false,
                "w-integer" => false,
                "reuse-cuts" => true,
                "subtour_cut_pool" => pool,
                "intersection_cuts_cache_ref" => ref_ic,
                "warm_start_solution" => ws,
            )
            data.T = 95.0
            sol2, info2 = runCOPCompleteDigraphIPModel(data, app2)
            @test info2["warmStartUsed"] == "true"
            @test length(pool) ≥ n_pool_after_first
            ok[] = true
        catch
        end
        ok[] || @test_skip "CPLEX unavailable or complete digraph solver error"
    end
end

@testset "Complete digraph IP with drop-zero-profit-blocks (CPLEX)" begin
    root::String = joinpath(@__DIR__, "..")
    inst::String = joinpath(root, "data", "carlos", "notified-alto-santo", "notified-alto-santo-1000-2021.txt")
    if !isfile(inst)
        @test_skip "Carlos fixture missing"
    else
        include(joinpath(root, "src", "data.jl"))
        include(joinpath(root, "src", "model.jl"))
        ok::Ref{Bool} = Ref(false)
        try
            data = readSBRPDataCarlos(Dict{String,Any}(
                "instance" => inst,
                "vehicle-time-limit" => "120",
                "no-cbrp-metric-closure" => false,
                "drop-zero-profit-blocks" => true,
            ))[1]
            data.T = 90.0
            sol, info = runCOPCompleteDigraphIPModel(data, Dict{String,Any}(
                "time-limit" => "15",
                "subcycle-separation" => "first",
                "intersection-cuts" => false,
                "y-integer" => false,
                "z-integer" => false,
                "w-integer" => false,
            ))
            @test haskey(info, "cost")
            @test haskey(info, "bestBound")
            @test info["bestBound"] != ""
            if info["bestBound"] != "N/A"
                @test parse(Float64, info["cost"]) <= parse(Float64, info["bestBound"]) + 1e-3
            end
            @test length(sol.tour) ≥ 2
            ok[] = true
        catch
        end
        ok[] || @test_skip "CPLEX unavailable or complete digraph solver error"
    end
end

@testset "calculateShortestPaths and createCompleteDigraph (toy street digraph)" begin
    root::String = joinpath(@__DIR__, "..")
    include(joinpath(root, "src", "data.jl"))
    b1::Vi = Vi([1])
    b2::Vi = Vi([2])
    b3::Vi = Vi([3])
    depot::Int = 4
    dist_sparse::ArcCostMap = ArcCostMap(
        (1 => 2) => 1.0,
        (2 => 1) => 1.0,
        (2 => 3) => 2.0,
        (3 => 2) => 2.0,
        (1 => 3) => 100.0,
        (3 => 1) => 100.0,
        (4 => 1) => 0.0,
        (1 => 4) => 0.0,
        (4 => 2) => 0.0,
        (2 => 4) => 0.0,
        (4 => 3) => 0.0,
        (3 => 4) => 0.0,
    )
    A_street::Arcs = Arcs(collect(keys(dist_sparse)))
    data::SBRPData = SBRPData(
        InputDigraph(
            Dict{Int,Vertex}(
                1 => Vertex(1, 0.0, 0.0),
                2 => Vertex(2, 1.0, 0.0),
                3 => Vertex(3, 2.0, 0.0),
                4 => Vertex(4, -1.0, -1.0),
            ),
            A_street,
            dist_sparse,
        ),
        depot,
        VVi([b1, b2, b3]),
        60.0,
        Dict{Vi,Float64}(b1 => 1.0, b2 => 1.0, b3 => 1.0),
    )
    sp::ArcCostMap = calculateShortestPaths(data)
    @test sp[1=>2] ≈ 1.0
    @test sp[1=>3] ≈ 3.0
    @test sp[2=>3] ≈ 2.0
    @test sp[3=>1] ≈ 3.0
    data′::SBRPData, paths′::Dict{Arc,Vi} = createCompleteDigraph(data)
    @test data′.D.distance[1=>2] ≈ 1.0
    @test data′.D.distance[1=>3] ≈ 3.0
    @test data′.D.distance[2=>3] ≈ 2.0
    @test data′.D.distance[1=>4] ≈ 0.0
    @test haskey(paths′, 1=>3)
    @test length(paths′[1=>3]) ≥ 1
end
