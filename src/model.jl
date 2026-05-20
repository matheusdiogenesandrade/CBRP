using CPLEX
using JuMP
using DataStructures
#using BenchmarkTools

const EPS = 1e-4

function unsetBinary(vars)
    unset_binary.(vars)
    set_lower_bound.(vars, zeros(length(vars)))
    set_upper_bound.(vars, ones(length(vars)))
end

function setBinary(vars)
    set_binary.(vars)
end

#=
Get all blocks with some node inside the component ´S´
input: 
- B::Vector{Vector{Int}} is the set of blocks
- S::Set{Int} is the component of nodes
=# 
function blocks(B::VVi, S::Si)::VVi
    return VVi(filter(block::Vi -> !isempty(intersect(S, block)), B))
end

include("ip_model.jl")
include("complete_subtour_cuts.jl")
include("path_cbrp_interfaces.jl")
include("brkga_path_warmstart.jl")
include("path_cbrp_ip.jl")
include("path_cbrp_subtour_cuts.jl")
