"""
    Algebra

Module implementing the algebra of hierarchical matrices. Includes the algebra between the various `AbstractMatrix` implemented in `HierarchicalMatrices` package.
"""
module Algebra

using ..HierarchicalMatrices.Matrices

import AbstractTrees
import LinearAlgebra

include("basics.jl")
include("svd.jl")

end
