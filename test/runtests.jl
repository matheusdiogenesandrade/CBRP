using Test

@testset "Carlos sparse reader vs metric closure" begin
    root::String = joinpath(@__DIR__, "..")
    inst::String = joinpath(root, "data", "carlos", "notified-alto-santo", "notified-alto-santo-1000-2021.txt")
    if !isfile(inst)
        @test_skip "Carlos fixture missing"
    else
        include(joinpath(root, "src", "data.jl"))
        app_dense = Dict{String, Any}(
            "instance" => inst,
            "vehicle-time-limit" => "120",
            "no-cbrp-metric-closure" => false,
        )
        app_sparse = Dict{String, Any}(
            "instance" => inst,
            "vehicle-time-limit" => "120",
            "no-cbrp-metric-closure" => true,
        )
        data_dense = readSBRPDataCarlos(app_dense)
        data_sparse = readSBRPDataCarlos(app_sparse)
        @test length(data_dense.D.A) > 0
        @test length(data_sparse.D.A) > 0
        @test length(data_dense.D.A) != length(data_sparse.D.A)
        @test length(data_dense.B) == length(data_sparse.B)
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
        sol::Union{Nothing, SBRPSolution} = nothing
        info::Union{Nothing, Dict{String, String}} = nothing
        try
            data = readSBRPDataCarlos(Dict{String, Any}(
                "instance" => inst,
                "vehicle-time-limit" => "120",
                "no-cbrp-metric-closure" => true,
            ))
            sol, info = runPathCbrpMipModel(data, Dict{String, Any}("time-limit" => "5"))
            ok[] = true
        catch
        end
        ok[] || @test_skip "CPLEX unavailable or Path-CBRP solver error"
        @test sol !== nothing && info !== nothing
        @test haskey(info, "cost")
        @test parse(Float64, info["cost"]) ≥ 0.0
        @test length(sol.tour) ≥ 2
    end
end
