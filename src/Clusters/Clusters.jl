"""
    Clusters

Module implementing various cluster trees and block cluster trees to be used in the construction of hierarchical matrices.
"""
module Clusters

import RecipesBase
import AbstractTrees
import Statistics

using ..HierarchicalMatrices: Maybe
using ..Geometry

export ClusterTree, BlockClusterTree

include("clustertree.jl")
include("splitter.jl")
include("blockclustertree.jl")
include("admissibility.jl")

end#module
