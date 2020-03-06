using SafeTestsets

@safetestset "Multiplication" begin
    using HierarchicalMatrices
    using HierarchicalMatrices: compress, ACA, PartialACA
    using HierarchicalMatrices.Geometry
    using HierarchicalMatrices.Clusters
    using LinearAlgebra
    @testset "mul!" begin
        let
            N    = 1000
            data = rand(Point{3,Float64},N)
            splitter   = Clusters.CardinalitySplitter(nmax=32)
            clt = ClusterTree(data,splitter;reorder=true)
            bclt = BlockClusterTree(clt,clt)
            comp = HierarchicalMatrices.ACA()
            f(x,y)::ComplexF64 = x==y ? 0.0 : exp(im*LinearAlgebra.norm(x-y))/LinearAlgebra.norm(x-y)
            M    = [f(x,y) for x in data, y in data]
            H    = HMatrix(M,bclt,comp)
            x    = rand(ComplexF64,N)
            @test norm(H*x-M*x) < comp.rtol*norm(M)
        end
    end
end
