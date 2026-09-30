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

"""
\\brief Minimum exact violation for a max-flow SEC to be reported (`--sec-min-violation`), shared by
the complete-digraph and Path-CBRP separations.

\\param app Parsed CLI options
\\return Threshold (default `1e-4`)
\\throws ArgumentError If the value is negative or not a number
"""
function secMinViolation(app::Dict{String,Any})::Float64
    raw = get(app, "sec-min-violation", 1e-4)
    v::Float64 = raw isa Real ? Float64(raw) : parse(Float64, String(raw))
    v >= 0.0 || throw(ArgumentError("sec-min-violation must be >= 0 (got $(v))"))
    return v
end

"""
\\brief Deepest B&B node at which fractional (user-cut) SECs are separated in the callback engine
(`--sec-user-cut-max-depth`), shared by the complete-digraph and Path-CBRP callbacks. Lazy SECs at
integer candidates are separated at every depth regardless.

\\param app Parsed CLI options
\\return Maximum depth (default `0` = root only; `-1` = unlimited)
\\throws ArgumentError If the value is below `-1` or not an integer
"""
function secUserCutMaxDepth(app::Dict{String,Any})::Int
    raw = get(app, "sec-user-cut-max-depth", 0)
    d::Int = raw isa Integer ? Int(raw) : parse(Int, String(raw))
    d >= -1 || throw(ArgumentError("sec-user-cut-max-depth must be >= -1 (got $(d))"))
    return d
end

"""
\\brief Whether fractional SECs may be separated at a node of the given depth.

\\param max_depth Limit from `secUserCutMaxDepth` (`-1` = unlimited)
\\param depth B&B depth of the current node (root = 0)
\\return `true` iff `max_depth < 0` or `depth <= max_depth`
"""
secUserCutDepthAllowed(max_depth::Int, depth::Int)::Bool = max_depth < 0 || depth <= max_depth

"""
\\brief B&B depth of the node a CPLEX generic callback was invoked at (root = 0).

\\param cb_data CPLEX generic callback context
\\return Node depth
\\throws ErrorException If CPLEX cannot report the depth
"""
function callbackNodeDepth(cb_data::CPLEX.CallbackContext)::Int
    depth_p = Ref{CPLEX.CPXINT}(0)
    status = CPLEX.CPXcallbackgetinfoint(cb_data, CPLEX.CPXCALLBACKINFO_NODEDEPTH, depth_p)
    status == 0 || error("CPXcallbackgetinfoint(NODEDEPTH) failed with status $(status)")
    return Int(depth_p[])
end

include("ip_model.jl")
include("complete_subtour_cuts.jl")
include("path_cbrp_interfaces.jl")
include("brkga_path_warmstart.jl")
include("path_cbrp_ip.jl")
include("path_cbrp_subtour_cuts.jl")
