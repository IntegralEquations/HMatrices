using HMatrices
using Test

using Random
Random.seed!(0)

include("rkmatrix_test.jl")
include("flexmatrix_test.jl")
include("rkflexmatrix_test.jl")
include("hmatrix_test.jl")
include("svd_test.jl")
include("compressor_test.jl")
include("blockdiagonal_test.jl")
include("preconditioner_test.jl")
# include("addition_test.jl")
# include("multiplication_test.jl")
# include("inverse_test.jl")
