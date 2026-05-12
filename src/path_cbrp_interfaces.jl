# Minimal accessors for Path-CBRP MILP on `SBRPData` (mirrors subset of external/CBRP interfaces.jl).

num_clusters(data::SBRPData)::Int = length(data.B)

get_cluster(data::SBRPData, idx::Int) = data.B[idx]

profit(data::SBRPData, idx::Int)::Float64 = data.profits[data.B[idx]]

service_time(data::SBRPData, idx::Int)::Float64 = blockTime(data, data.B[idx])

arc_time(data::SBRPData, arc::Arc)::Float64 = time(data, arc)

depot_node(data::SBRPData)::Int = data.depot
