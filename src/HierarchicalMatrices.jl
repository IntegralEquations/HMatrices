module HierarchicalMatrices

import LinearAlgebra

export HMatrix

include("Interfaces.jl")
include("utils.jl")
include("Geometry/Geometry.jl")
include("Clusters/Clusters.jl")
include("Matrices/Matrices.jl")
include("Algebra/Algebra.jl")

end # module
