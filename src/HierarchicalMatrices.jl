module HierarchicalMatrices

import LinearAlgebra
import AbstractTrees

export HMatrix, LazyMatrix, Geometry, Clusters

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
