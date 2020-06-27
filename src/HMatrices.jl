module HMatrices

import LinearAlgebra: rank, mul!, svd, svd!, norm, axpby!, axpy!, rmul!, inv!, lu, lu!, ldiv!, rdiv!, triu, tril
import Base: +, -, *, inv, copy
import AbstractTrees: children
import SparseArrays: sparse, SparseMatrixCSC

using LinearAlgebra: Adjoint, UnitLowerTriangular, LowerTriangular, UnitUpperTriangular, UpperTriangular,
                     I, SVD, qr, qr!, adjoint!, LU, full!, Diagonal, dot, UniformScaling, det, diag
using SparseArrays: SparseMatrixCSC, AbstractSparseArray, rowvals, nonzeros, nzrange
using Base.Threads: @spawn
using ComputationalResources: CPU1, CPUThreads, AbstractResource
using AbstractTrees: TreeIterator, Leaves, PreOrderDFS
using RecipesBase: @series, @recipe
using RecipesBase
using UnsafeArrays: uview

export
    # Types
    HMatrix,
    LazyMatrix,
    CPUThreads

# Constants
const Maybe{T}  = Union{Tuple{},T}

# Includes
include("interface.jl")
include("utils.jl")
include("lazymatrix.jl")
include("flexmatrix.jl")
include("rkmatrix.jl")
include("hmatrix.jl")
include("conversion.jl")
include("adjoint.jl")
include("svd.jl")
include("norm.jl")
include("compressor.jl")
include("addition.jl")
include("multiplication.jl")
include("inverse.jl")
include("triangular.jl")
include("lu.jl")

include("parameters.jl")

end # module
