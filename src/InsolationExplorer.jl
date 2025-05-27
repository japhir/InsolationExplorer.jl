module InsolationExplorer

using FileIO
using CSV
using DataFrames
using Downloads: download
using ColorTypes
using GLMakie
using GeometryBasics
using Rotations

include("insolation.jl")
include("explore_insolation.jl")
include("explore_solution.jl")

export insolation
export explore_insolation
export explore_solution

end # module
