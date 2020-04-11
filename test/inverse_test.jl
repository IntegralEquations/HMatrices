using SafeTestsets

@safetestset "Inverse" begin
    using HierarchicalMatrices
    using Clusters
    using HierarchicalMatrices: PartialACA, RkMatrix
    using LinearAlgebra

    N,r    = 500,4
    T    = ComplexF64
    data = rand(Clusters.Point{2,Float64},N)
    splitter   = Clusters.GeometricMinimalSplitter()
    clt  = Clusters.ClusterTree(data,splitter;reorder=true)
    adm  = Clusters.StrongAdmissibilityStd(100)
    bclt = Clusters.BlockTree(clt,clt,adm)
    comp = HierarchicalMatrices.PartialACA(rtol=1e-6)
    f(x,y)::ComplexF64 = exp(-norm(x-y))
    L    = LazyMatrix(f,data,data)
    H    = HMatrix(L,bclt,comp)
    _H   = Matrix(H)
    R    = RkMatrix(rand(T,N,r),rand(T,N,r))
    _R   = Matrix(R)
    M    = rand(T,N,N)
    _M   = Matrix(M)
    @testset "inv" begin
        #TODO: find a better test for the inverse. The one below is poor due to the
        #condition number of H
        Hinv = inv(H) |> Matrix
        @test norm(Hinv - inv(_H)) < 1e-4
    end
    @testset "ldiv!" begin
        X    = rand(T,N,5)
        L  = UnitLowerTriangular(H)
        _L = UnitLowerTriangular(_H)
        @test ldiv!(L,copy(X)) ≈ _L\X
        @test ldiv!(L,deepcopy(R)) ≈ _L\_R
        @test ldiv!(L,deepcopy(H)) ≈ _L\_H
        L  = LowerTriangular(H)
        _L = LowerTriangular(_H)
        @test ldiv!(L,copy(X)) ≈ _L\X
        @test ldiv!(L,deepcopy(R)) ≈ _L\_R
        @test ldiv!(L,deepcopy(H)) ≈ _L\_H
        U  = UnitUpperTriangular(H)
        _U = UnitUpperTriangular(_H)
        @test ldiv!(U,copy(X)) ≈ _U\X
        @test ldiv!(U,deepcopy(R)) ≈ _U\_R
        @test ldiv!(U,deepcopy(H)) ≈ _U\_H
        U  = UpperTriangular(H)
        _U = UpperTriangular(_H)
        @test ldiv!(U,copy(X)) ≈ _U\X
        @test ldiv!(U,deepcopy(R)) ≈ _U\_R
        @test ldiv!(U,deepcopy(H)) ≈ _U\_H
    end
    @testset "rdiv!" begin
        X    = rand(T,5,N)
        L  = UnitLowerTriangular(H)
        _L = UnitLowerTriangular(_H)
        @test rdiv!(copy(X),L) ≈ X/_L
        @test rdiv!(deepcopy(R),L) ≈ _R/_L
        @test rdiv!(deepcopy(H),L) ≈ _H/_L
        L  = LowerTriangular(H)
        _L = LowerTriangular(_H)
        @test rdiv!(copy(X),L) ≈ X/_L
        @test rdiv!(deepcopy(R),L) ≈ _R/_L
        @test rdiv!(deepcopy(H),L) ≈ _H/_L
        U  = UnitUpperTriangular(H)
        _U = UnitUpperTriangular(_H)
        @test rdiv!(copy(X),U) ≈ X/_U
        @test rdiv!(deepcopy(R),U) ≈ _R/_U
        @test rdiv!(deepcopy(H),U) ≈ _H/_U
    end
end
