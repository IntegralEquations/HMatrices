module HierarchicalMatrices

using Base.Threads: @spawn
using LinearAlgebra
using ComputationalResources
using AbstractTrees
using RecipesBase
using UnsafeArrays

import AbstractTrees: children

import LinearAlgebra: rank, mul!, svd, svd!, norm, axpby!, axpy!, rmul!, inv!, lu, lu!, ldiv!, rdiv!
import Base: +, -, *, inv

export HMatrix, LazyMatrix, Clusters
export CPUThreads

include("Interface.jl")
include("Parameters.jl")
include("utils.jl")
include("Clusters/Clusters.jl")

using .Clusters
using .Parameters

import .Interface: rowrange, colrange, getchildren, getparent, isleaf, isroot, isadmissible

const Maybe{T}  = Union{Tuple{},T}

include("Matrices/Matrices.jl")
include("Algebra/Algebra.jl")

end # module
