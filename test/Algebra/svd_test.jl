using SafeTestsets

@safetestset "RkMatrix" begin
    using HierarchicalMatrices.Matrices
    using LinearAlgebra: svd, svd!
    m,n,r = 20, 10, 5
    Rk = rand(RkMatrix,m,n,r)
    F = svd(Rk)
    tmp = svd(Matrix(Rk))
    @test Matrix(F) ≈ Matrix(tmp)
    @test F.S[1:r] ≈ tmp.S[1:r]
    F = svd!(copy(Rk)) # changes Rk
    @test Matrix(F) ≈ Matrix(tmp)
    @test F.S[1:r] ≈ tmp.S[1:r]
end
