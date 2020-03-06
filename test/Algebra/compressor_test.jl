using SafeTestsets

@safetestset "Compressor" begin
    using HierarchicalMatrices
    using HierarchicalMatrices: ACA, PartialACA,compress
    using LinearAlgebra: norm, Diagonal, rank
    T = ComplexF64
    m,n,r = 100,100,100
    S = Diagonal(exp.(-(1:r)))
    M  = 1e3.*rand(T,m,r)*S*rand(T,r,n)
    irange,jrange = 1:size(M,1),1:size(M,2)
    @testset "aca_full" begin
        atol = 1e-5
        aca = ACA(atol=atol)
        R = compress(M,irange,jrange,aca)
        @info  rank(R)
        norm(Matrix(R) - M)
        @test norm(Matrix(R) - M) < atol
        rtol = 1e-5
        aca = ACA(rtol=rtol)
        R = compress(M,irange,jrange,aca)
        norm(Matrix(R) - M)
        @test norm(Matrix(R) - M) < rtol*norm(M)
        r = 10
        aca = ACA(rank=r)
        R = compress(M,irange,jrange,aca)
        @test rank(R) == r
    end
    @testset "aca_partial" begin
        atol = 1e-5
        aca = PartialACA(atol=atol)
        R = compress(M,irange,jrange,aca)
        @info  rank(R)
        @test norm(Matrix(R) - M) < atol
        rtol = 1e-5
        aca = ACA(rtol=rtol)
        R = compress(M,irange,jrange,aca)
        @test norm(Matrix(R) - M) < rtol*norm(M)
        r = 10
        aca = ACA(rank=r)
        R = compress(M,irange,jrange,aca)
        @test rank(R) == r
    end
end
