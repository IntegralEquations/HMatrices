using SafeTestsets

@safetestset "ClusterTree" begin
    using HierarchicalMatrices.Clusters
    using HierarchicalMatrices.Clusters: isroot, isleaf, Point, bounding_box
    @testset "2d binary clusters" begin
        let
            data = rand(Point{2,Float64},1000)
            splitter   = Clusters.GeometricSplitter(nmax=32)
            clt = ClusterTree(data,splitter;reorder=false)
            @test clt.data == data[clt.perm]
            splitter   = Clusters.GeometricMinimalSplitter(nmax=32)
            clt = ClusterTree(data,splitter;reorder=false)
            @test clt.data == data[clt.perm]
            splitter   = Clusters.CardinalitySplitter(nmax=32)
            clt = ClusterTree(data,splitter;reorder=false)
            @test clt.data == data[clt.perm]
        end
    end

    @testset "3d" begin
        let
            data = rand(Point{3,Float64},1000)
            splitter   = Clusters.GeometricSplitter(nmax=32)
            clt = ClusterTree(data,splitter;reorder=false)
            @test clt.data == data[clt.perm]
            splitter   = Clusters.GeometricMinimalSplitter(nmax=32)
            clt = ClusterTree(data,splitter;reorder=false)
            @test clt.data == data[clt.perm]
            splitter   = Clusters.CardinalitySplitter(nmax=32)
            clt = ClusterTree(data,splitter;reorder=false)
            @test clt.data == data[clt.perm]
        end
    end
end
