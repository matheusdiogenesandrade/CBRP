include("symbols.jl")

import Unicode

mutable struct Vertex
   id_vertex::Int
   pos_x::Float64
   pos_y::Float64
end

mutable struct InputDigraph # Directed graph
   V::Dict{Int, Vertex} # set of vertices
   A::Arcs              # set of arcs
   distance::ArcCostMap # time in 40 KM/h
end

δ⁺(A::Arcs, i::Int)::Arcs = filter((j, k)::Arc -> j == i, A)
δ⁺(A::Arcs, S::Si)::Arcs  = filter((i, j)::Arc -> i in S && !in(j, S), A)
δ⁺(A::Arcs, S::Vi)::Arcs  = filter((i, j)::Arc -> i in S && !in(j, S), A)

δ⁺(A::ArcsSet, i::Int)::ArcsSet = filter((j, k)::Arc -> j == i, A)
δ⁺(A::ArcsSet, S::Si)::ArcsSet  = filter((i, j)::Arc -> i in S && !in(j, S), A)
δ⁺(A::ArcsSet, S::Vi)::ArcsSet  = filter((i, j)::Arc -> i in S && !in(j, S), A)

δ⁻(A::Arcs, i::Int)::Arcs = filter((j, k)::Arc -> k == i, A)
δ⁻(A::Arcs, S::Si)::Arcs  = filter((i, j)::Arc -> j in S && !in(i, S), A)
δ⁻(A::Arcs, S::Vi)::Arcs  = filter((i, j)::Arc -> j in S && !in(i, S), A)

δ⁻(A::ArcsSet, i::Int)::ArcsSet = filter((j, k)::Arc -> k == i, A)
δ⁻(A::ArcsSet, S::Si)::ArcsSet  = filter((i, j)::Arc -> j in S && !in(i, S), A)
δ⁻(A::ArcsSet, S::Vi)::ArcsSet  = filter((i, j)::Arc -> j in S && !in(i, S), A)

d⁺(A::ArcsSet, i::Int)::Int64 = length(δ⁺(A, i))
d⁻(A::ArcsSet, i::Int)::Int64 = length(δ⁻(A, i))
d⁺(A::Arcs, i::Int)::Int64    = length(δ⁺(A, i))
d⁻(A::Arcs, i::Int)::Int64    = length(δ⁻(A, i))

v⁺(A::ArcsSet, i::Int)::Vi = map((i, j)::Arc -> j, δ⁺(A, i))
v⁻(A::ArcsSet, i::Int)::Vi = map((i, j)::Arc -> i, δ⁻(A, i))

v⁺(A::Arcs, i::Int)::Vi = map((i, j)::Arc -> j, δ⁺(A, i))
v⁻(A::Arcs, i::Int)::Vi = map((i, j)::Arc -> i, δ⁻(A, i))

include("sbrp.jl")
include("simplify_street_graph.jl")
