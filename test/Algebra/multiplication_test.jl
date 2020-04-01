@safetestset "Multiplication" begin
    using HierarchicalMatrices
    using HierarchicalMatrices: ACA, PartialACA
    using LinearAlgebra

    N    = 1000
    data = rand(Clusters.Point{2,Float64},N)
    splitter   = Clusters.GeometricMinimalSplitter()
    clt  = Clusters.ClusterTree(data,splitter;reorder=true)
    bclt = Clusters.BlockTree(clt,clt)
    comp = HierarchicalMatrices.PartialACA(rtol=1e-6)
    f(x,y)::ComplexF64 = x==y ? 0.0 : exp(im*LinearAlgebra.norm(x-y))/LinearAlgebra.norm(x-y)
    M    = LazyMatrix(f,data,data)
    H    = HMatrix(M,bclt,comp)
    _H   = Matrix(H)
    a    = rand(); b = rand();
    @testset "CPU1 mul!" begin
        x    = rand(ComplexF64,N)
        y    = similar(x)
        @test mul!(y,H,x) ≈ _H*x
        @test H*x ≈ _H*x
        tmp = b*y + a*_H*x
        @test mul!(y,H,x,a,b) ≈ tmp
    end
    @testset "CPUThreads mul!" begin
        x    = rand(ComplexF64,N)
        y    = zero(x)
        @test mul!(CPUThreads(),y,H,x,true,false) ≈ _H*x
        tmp = b*y + a*_H*x
        @test mul!(CPUThreads(),y,H,x,a,b) ≈ tmp
    end
end
