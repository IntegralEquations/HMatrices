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

export FlexMatrix, RkFlexMatrix, RkMatrix, HMatrix

include("rkmatrix.jl")
include("flexmatrix.jl")
include("rkflexmatrix.jl")


end#module
