using SafeTestsets

@safetestset "HMatrix" begin
    using HierarchicalMatrices
    using HierarchicalMatrices: compress, ACA, PartialACA
    using HierarchicalMatrices.Geometry
    using HierarchicalMatrices.Clusters
    using LinearAlgebra
    @testset "Assembly" begin
        let
            N    = 100
            data = rand(Point{3,Float64},N)
            splitter   = Clusters.CardinalitySplitter(nmax=32)
            clt = ClusterTree(data,splitter;reorder=true)
            bclt = BlockClusterTree(clt,clt)
            comp = HierarchicalMatrices.ACA()
            f(x,y)::ComplexF64 = x==y ? 0.0 : exp(im*LinearAlgebra.norm(x-y))/LinearAlgebra.norm(x-y)
            M    = [f(x,y) for x in data, y in data]
            H    = HMatrix(M,bclt,comp)
            @test norm(H-M,2) < comp.rtol*norm(M)
        end
    end
end
