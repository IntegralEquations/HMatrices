"""
    Matrices

Module implementing various matrices types involved in the  construction of Hierarchical Matrices.
"""
module Matrices

import RecipesBase
import AbstractTrees
import LinearAlgebra

using ..HierarchicalMatrices: Maybe
using ..HierarchicalMatrices.Clusters

export RkMatrix, HMatrix

include("rkmatrix.jl")

end#module
