using SafeTestsets

@safetestset "Inverse" begin
    using HierarchicalMatrices
    using HierarchicalMatrices: PartialACA, RkMatrix
    using LinearAlgebra

    #TODO: find a better test for the inverse. The one below is poor due to the
    #condition number of H
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
    X    = rand(T,N,5)
    @testset "inv" begin
        Hinv = inv(H) |> Matrix
        @test norm(Hinv - inv(_H)) < 1e-3
    end
    @testset "ldiv!" begin
        L  = UnitLowerTriangular(H)
        _L = UnitLowerTriangular(_H)
        tmp = copy(X)
        ldiv!(L,tmp)
        @test tmp ≈ inv(_L)*X
        tmp = copy(R)
        ldiv!(L,tmp)
        @test tmp ≈ inv(_L)*_R
        tmp = copy(H)
        ldiv!(L,tmp)
        @test tmp ≈ inv(_L)*_H
        U  = UpperTriangular(H)
        _U = UpperTriangular(_H)
        tmp = copy(X)
        ldiv!(U,tmp)
        @test tmp ≈ inv(_U)*X        
        tmp = copy(R)
        ldiv!(U,tmp)
        @test tmp ≈ inv(_U)*R
        tmp = copy(H)
        ldiv!(U,tmp)
        @test tmp ≈ inv(_U)*_H
    end
    @testset "rdiv!" begin
        L  = UnitLowerTriangular(H)
        _L = UnitLowerTriangular(_H)
        tmp = copy(X)
        ldiv!(L,tmp)
        @test tmp ≈ inv(_L)*X
        tmp = copy(R)
        ldiv!(L,tmp)
        @test tmp ≈ inv(_L)*_R
        tmp = copy(H)
        ldiv!(L,tmp)
        @test tmp ≈ inv(_L)*_H
        U  = UpperTriangular(H)
        _U = UpperTriangular(_H)
        tmp = copy(X)
        ldiv!(U,tmp)
        @test tmp ≈ inv(_U)*X        
        tmp = copy(R)
        ldiv!(U,tmp)
        @test tmp ≈ inv(_U)*R
        tmp = copy(H)
        ldiv!(U,tmp)
        @test tmp ≈ inv(_U)*_H
    end
end
