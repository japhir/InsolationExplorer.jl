module InsolationExplorer

using FileIO
using Downloads: download
using ColorTypes
using GLMakie
using GeometryBasics
using Rotations

include("insolation.jl")
include("explore_insolation.jl")
# include("explore_solution.jl") # TODO

export insolation
export explore_insolation
# export explore_solution # TODO

end # module
