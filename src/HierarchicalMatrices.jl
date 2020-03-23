module HierarchicalMatrices

using Base.Threads: @spawn
using LinearAlgebra
using ComputationalResources
using AbstractTrees

import AbstractTrees: children

export HMatrix, LazyMatrix, Clusters
export CPUThreads

include("Interface.jl")
include("Parameters.jl")
include("utils.jl")
include("Clusters/Clusters.jl")

using .Clusters
using .Parameters

import .Interface: rowrange, colrange, getchildren, getparent, isleaf, isroot, isadmissible

include("Matrices/Matrices.jl")
include("Algebra/Algebra.jl")

include("treeutils.jl")

end # module
