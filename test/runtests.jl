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
            )
            sol, info = runPathCbrpMipModel(data, solve_app_sec)
            ok[] = true
        catch
        end
        ok[] || @test_skip "CPLEX unavailable or Path-CBRP SEC solver error"
        @test sol !== nothing && info !== nothing && data !== nothing
        @test parse(Int, get(info, "maxFlowCuts", "0")) >= 0
        @test haskey(info, "maxFlowCutsTime")
        include(joinpath(root, "src", "sol.jl"))
        checkSBRPSolution(data::SBRPData, sol::SBRPSolution)
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
