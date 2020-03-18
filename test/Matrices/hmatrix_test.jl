using SafeTestsets

@safetestset "HMatrix" begin
    using HierarchicalMatrices
    using HierarchicalMatrices: ACA, PartialACA
    using HierarchicalMatrices.Geometry
    using HierarchicalMatrices.Clusters
    using LinearAlgebra
    @testset "Assembly" begin
        let
            N    = 1000
            data = rand(Point{2,Float64},N)
            splitter   = Clusters.CardinalitySplitter(nmax=32)
            clt = ClusterTree(data,splitter;reorder=true)
            bclt = BlockClusterTree(clt,clt)
            f(x,y)::ComplexF64 = x==y ? 0.0 : exp(im*LinearAlgebra.norm(x-y))/LinearAlgebra.norm(x-y)
            M    = LazyMatrix(f,data,data)
            comp = HierarchicalMatrices.ACA()
            H    = HMatrix(M,bclt,comp)
            @test norm(H-M,2) < comp.rtol*norm(M)
            comp = HierarchicalMatrices.PartialACA()
            H    = HMatrix(M,bclt,comp)
            @test norm(H-M,2) < comp.rtol*norm(M)
        end
    end
end
