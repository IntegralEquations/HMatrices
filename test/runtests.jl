using HierarchicalMatrices
using Test

using Random
Random.seed!(0)

include("Clusters/runtests.jl")
include("Matrices/runtests.jl")
include("Algebra/runtests.jl")
