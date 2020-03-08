"""
    Clusters

Module implementing various cluster trees and block cluster trees to be used in the construction of hierarchical matrices.
"""
module Clusters

import RecipesBase
import AbstractTrees
import Statistics

using  ..Parameters

import ..Interfaces: container, split, diameter, distance, rowrange, colrange,
                     isleaf, isroot, getchildren, getparent, isadmissible

using  ..HierarchicalMatrices: Maybe

export ClusterTree, BlockClusterTree

include("clustertree.jl")
include("splitter.jl")
include("blockclustertree.jl")
include("admissibility.jl")

end#module
