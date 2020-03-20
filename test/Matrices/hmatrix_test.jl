using SafeTestsets

@safetestset "HMatrix" begin
    using HierarchicalMatrices
    using HierarchicalMatrices: ACA, PartialACA
    using ComputationalResources
    using LinearAlgebra
    N    = 1000
    data = rand(Geometry.Point{2,Float64},N)
    splitter   = Clusters.CardinalitySplitter(nmax=128)
    clt  = Clusters.ClusterTree(data,splitter;reorder=true)
    bclt = Clusters.BlockClusterTree(clt,clt)
    f(x,y)::ComplexF64 = x==y ? 0.0 : exp(im*LinearAlgebra.norm(x-y))/LinearAlgebra.norm(x-y)
    M    = LazyMatrix(f,data,data)
    comp = HierarchicalMatrices.ACA(rtol=1e-6)
    @testset "Assembly CPU1" begin
            H    = HMatrix(M,bclt,comp)
            @test norm(H-M,2) < comp.rtol*norm(M)
            comp = HierarchicalMatrices.PartialACA()
            H    = HMatrix(M,bclt,comp)
            @test norm(H-M,2) < comp.rtol*norm(M)
    end
    @testset "Assembly CPUThreads" begin
        H    = HMatrix(CPUThreads(),M,bclt,comp)
        @test norm(H-M,2) < comp.rtol*norm(M)
        comp = HierarchicalMatrices.PartialACA()
        H    = HMatrix(M,bclt,comp)
        @test norm(H-M,2) < comp.rtol*norm(M)
    end
end
