using SafeTestsets

@safetestset "Compressor" begin
    using HierarchicalMatrices
    using HierarchicalMatrices: compress, ACA, PartialACA, TSVD
    using LinearAlgebra: norm, Diagonal, rank
    T = ComplexF64
    m,n,r = 100,100,100
    S = Diagonal(exp.(-(1:r)))
    M  = 1e3.*rand(T,m,r)*S*rand(T,r,n)
    irange,jrange = 1:size(M,1),1:size(M,2)
    @testset "aca_full" begin
        atol = 1e-5
        aca = ACA(atol=atol)
        R = aca(M,irange,jrange)
        norm(Matrix(R) - M)
        @test norm(Matrix(R) - M) < atol
        rtol = 1e-5
        aca = ACA(rtol=rtol)
        R = aca(M,irange,jrange)
        norm(Matrix(R) - M)
        @test norm(Matrix(R) - M) < rtol*norm(M)
        r = 10
        aca = ACA(rank=r)
        R = aca(M,irange,jrange)
        @test rank(R) == r
    end
    @testset "aca_partial" begin
        atol = 1e-5
        aca = PartialACA(atol=atol)
        R = aca(M,irange,jrange)
        @test norm(Matrix(R) - M) < atol
        rtol = 1e-5
        aca = PartialACA(rtol=rtol)
        R = aca(M,irange,jrange)
        @test norm(Matrix(R) - M) < rtol*norm(M)
        r = 10
        aca = PartialACA(rank=r)
        R = aca(M,irange,jrange)
        @test rank(R) == r
        using HierarchicalMatrices: FlexMatrix, RkFlexMatrix, pushcross!, _update_frob_norm
        m,n,r = 10,10,4
        T = ComplexF64
        # T = Float64
        A = rand(T,m,r)
        B = rand(T,n,r)
        R = RkFlexMatrix(A,B)
        old_norm = norm(R)
        a = rand(T,m)
        b = rand(T,n)
        pushcross!(R,a,b)
        new_norm = norm(R)
        @test new_norm ≈ _update_frob_norm(old_norm,R)
    end
    @testset "truncated svd" begin
        # atol  = 1e-5
        # tsvd  = TSVD(atol=atol)
        # tmp   = compress(R,tsvd)
        # @test norm(Matrix(R) - M) < atol
        # rtol = 1e-5
        # aca = PartialACA(rtol=rtol)
        # R = aca(M,irange,jrange)
        # @test norm(Matrix(R) - M) < rtol*norm(M)
        # r = 10
        # aca = PartialACA(rank=r)
        # R = aca(M,irange,jrange)
        # @test rank(R) == r
    end
end
