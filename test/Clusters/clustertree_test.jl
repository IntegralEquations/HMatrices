using SafeTestsets

using HierarchicalMatrices
using Test

@safetestset "ClusterTree" begin
    using HierarchicalMatrices.Clusters
    using HierarchicalMatrices.Clusters: isroot, isleaf
    using HierarchicalMatrices.Geometry: Point, bounding_box
    @testset "2d" begin
        let
            data = rand(Point{2,Float64},1000)
            split_type = Clusters.GeometricSplit
            n_el       = 100
            reorder    = false
            clt_opts   = Clusters.ClusterTreeOptions(split_type,n_el,reorder)
            clt = ClusterTree(data,clt_opts)
            @test clt.data == data[clt.perm]
            Clusters.diameter(clt)
        end
    end

    @testset "3d" begin
        let
            data = rand(Point{3,Float64},1000)
            split_type = Clusters.GeometricMinimalSplit
            n_el       = 100
            reorder    = false
            clt_opts = Clusters.ClusterTreeOptions(split_type,n_el,reorder)
            clt = ClusterTree(data,clt_opts)
            @test clt.data == data[clt.perm]
        end
    end
end
