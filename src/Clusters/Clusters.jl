"""
    Clusters

Module implementing various cluster trees and block cluster trees to be used in the construction of hierarchical matrices.
"""
module Clusters

using LinearAlgebra: norm
using Statistics: median
using AbstractTrees
using RecipesBase

import AbstractTrees: children

using  ..Parameters

import ..Interface: rowrange, colrange, isleaf, isroot, getchildren, getparent, isadmissible

export ClusterTree, BlockTree, Point

const Maybe{T} = Union{Tuple{},T}

include("point.jl")
include("hyperrectangle.jl")
include("clustertree.jl")
include("splitter.jl")
include("blocktree.jl")
include("admissibility.jl")

end#module
