"""
    Clusters

Module implementing various cluster trees and block cluster trees to be used in the construction of hierarchical matrices.
"""
module Clusters

import RecipesBase

using ..HierarchicalMatrices: Maybe
using ..Geometry

export ClusterTree, BlockClusterTree

include("clustertree.jl")
# include("blockclustertree.jl")
# include("admissibilityfunction.jl")

end#module
