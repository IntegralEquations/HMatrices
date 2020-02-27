module HierarchicalMatrices

################################################################################
## GLOBAL CONSTANTS
################################################################################
const Maybe{T} = Union{Tuple{},T}

include("Geometry/Geometry.jl")
include("Clusters/Clusters.jl")
include("Matrices/Matrices.jl")
include("Algebra/Algebra.jl")

end # module
