module HierarchicalMatrices

using Base.Threads: @spawn
using LinearAlgebra
using AbstractTrees

import LinearAlgebra
import AbstractTrees

using ComputationalResources

export HMatrix, LazyMatrix, Geometry, Clusters
export CPUThreads

include("Interfaces.jl")
include("Parameters.jl")
include("utils.jl")
include("Geometry/Geometry.jl")
include("Clusters/Clusters.jl")

using .Clusters
using .Parameters

import .Interfaces: rowrange, colrange, getchildren, getparent, isleaf, isroot, isadmissible

include("Matrices/Matrices.jl")
include("Algebra/Algebra.jl")

end # module
