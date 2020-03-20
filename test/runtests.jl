using HierarchicalMatrices
using Test

using Random
Random.seed!(0)

include("Geometry/runtests.jl")
include("Clusters/runtests.jl")
include("Matrices/runtests.jl")
include("Algebra/runtests.jl")
