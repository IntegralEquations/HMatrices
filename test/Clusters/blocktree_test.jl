using SafeTestsets

@safetestset "BlockTree" begin
    using HierarchicalMatrices.Clusters
    using HierarchicalMatrices.Clusters: isadmissible, AdmissibilityStandard, Point
    let
        Xdata    = rand(Point{2,Float64},1000)
        adm_fun  = AdmissibilityStandard(2)
        Xclt     = ClusterTree(Xdata;reorder=false)
        bclt = BlockTree(Xclt,Xclt,adm_fun)
        @test isadmissible(bclt) == false
        adm_fun = AdmissibilityStandard(4)
        Ydata = Point{2,Float64}.(rand(10000).+10,rand(10000).+10)
        Yclt  = ClusterTree(Ydata;reorder=false)
        bclt = BlockTree(Xclt,Yclt,adm_fun)
        @test isadmissible(bclt) == true
    end
end
