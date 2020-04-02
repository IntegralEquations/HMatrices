using SafeTestsets

@safetestset "Inverse" begin
    using HierarchicalMatrices
    using HierarchicalMatrices: PartialACA
    using LinearAlgebra

    #TODO: find a better test for the inverse. The one below is poor due to the
    #condition number of H
    N    = 500
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
    Hinv = inv(H) |> Matrix
    @test norm(Hinv - inv(_H)) < 1e-3
end
