"""
    Clusters

Module implementing various cluster trees and block cluster trees to be used in the construction of hierarchical matrices.
"""
module Clusters

using LinearAlgebra: norm
using RecipesBase: @recipe, @series
import AbstractTrees
import Statistics

using  ..Parameters

import ..Interface: rowrange, colrange, isleaf, isroot, getchildren, getparent, isadmissible

using  ..HierarchicalMatrices: Maybe

export ClusterTree, BlockTree, Point

include("point.jl")
include("utils.jl")
include("hyperrectangle.jl")
include("clustertree.jl")
include("splitter.jl")
include("blocktree.jl")
include("admissibility.jl")

end#module
