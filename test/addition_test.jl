using SafeTestsets

@safetestset "Addition" begin
    using HierarchicalMatrices
    using Clusters
    using Clusters: Point
    using HierarchicalMatrices: ACA, PartialACA, RkMatrix
    using LinearAlgebra

    N,r     = 500,3
    atol    = 1e-6
    data    = rand(Point{3,Float64},N)
    splitter = GeometricMinimalSplitter()
    clt   = ClusterTree(data,splitter;reorder=true)
    adm  = Clusters.StrongAdmissibilityStd(Inf)
    bclt = Clusters.BlockTree(clt,clt,adm)
    comp = HierarchicalMatrices.PartialACA(atol=atol)
    f(x,y)::ComplexF64 = x==y ? 0.0 : exp(im*LinearAlgebra.norm(x-y))/LinearAlgebra.norm(x-y)
    M    = LazyMatrix(f,data,data)
    T    = ComplexF64
    H    = HMatrix(M,bclt,comp)
    _H   = Matrix(H)
    R    = RkMatrix(rand(T,N,r),rand(T,N,r))
    _R   = Matrix(R)
    M    = rand(T,N,N)
    _M   = Matrix(M)
    a    = rand()
    b    = rand()
    @testset "+" begin
        @test M+R ≈ _M + _R
        @test M+H ≈ _M + _H
        @test R+M ≈ _R + _M
        @test R+R ≈ _R + _R
        @test R+H ≈ _R + _H
        @test H+M ≈ _H + _M
        @test H+R ≈ _H + _R
        @test H+H ≈ _H + _H
    end
    @testset "axpby!" begin
        @test axpby!(a,M,b,deepcopy(R)) ≈ a*_M + b*_R
        @test axpby!(a,M,b,deepcopy(H)) ≈ a*_M + b*_H
        @test axpby!(a,R,b,deepcopy(M)) ≈ a*_R + b*_M
        @test axpby!(a,R,b,deepcopy(R)) ≈ a*_R + b*_R
        @test axpby!(a,R,b,deepcopy(H)) ≈ a*_R + b*_H
        @test axpby!(a,H,b,deepcopy(M)) ≈ a*_H + b*_M
        @test axpby!(a,H,b,deepcopy(R)) ≈ a*_H + b*_R
        @test axpby!(a,H,b,deepcopy(H)) ≈ a*_H + b*_H
    end
end
