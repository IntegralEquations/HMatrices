"""
    Clusters

Module implementing various cluster trees and block cluster trees to be used in the construction of hierarchical matrices.
"""
module Clusters

import RecipesBase
import AbstractTrees
import Statistics

import ..Interfaces: container, split, diameter, distance
using  ..HierarchicalMatrices: Maybe

export ClusterTree, BlockClusterTree

include("clustertree.jl")
include("splitter.jl")
include("blockclustertree.jl")
include("admissibility.jl")

end#module
