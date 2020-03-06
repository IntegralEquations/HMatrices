using SafeTestsets

@safetestset "aca" begin
    using HierarchicalMatrices.Matrices
    using HierarchicalMatrices.Algebra: aca_full, aca_partial
    using LinearAlgebra: norm, Diagonal
    @testset "aca_full" begin
        m,n,r = 100,100,30
        S = Diagonal(exp.(-(1:r)))
        M  = 1e5.*rand(m,r)*S*rand(r,n)
        atol = 1e-5
        F  = aca_full(M;atol=atol)
        @time aca_full(M;atol=atol);
        @test norm(Matrix(F) - M) < atol
        rtol = 1e-5
        F  = aca_full(M;rtol=rtol)
        @test !(norm(Matrix(F) - M) < rtol)
        @test norm(Matrix(F) - M)   < norm(M)*rtol
    end
    @testset "aca_partial" begin
        m,n,r = 100,100,30
        S = Diagonal(exp.(-(1:r)))
        M  = 1e5.*rand(m,r)*S*rand(r,n)
        atol = 1e-5
        F  = aca_partial(M;atol=atol)
        @time aca_partial(M;atol=atol);
        @test norm(Matrix(F) - M) < atol
        rtol = 1e-5
        F  = aca_full(M;rtol=rtol)
        @test !(norm(Matrix(F) - M) < rtol)
        @test norm(Matrix(F) - M)   < norm(M)*rtol
    end
end
