module HierarchicalMatrices

################################################################################
## GLOBAL CONSTANTS
################################################################################
const Maybe{T} = Union{Tuple{},T}

function debug(flag=true)
    if flag
        @eval ENV["JULIA_DEBUG"] = "HierarchicalMatrices"
    else
        @eval ENV["JULIA_DEBUG"] = ""
    end
end

include("Geometry/Geometry.jl")
include("Clusters/Clusters.jl")
include("Matrices/Matrices.jl")
include("Algebra/Algebra.jl")

end # module
