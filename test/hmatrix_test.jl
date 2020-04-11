using SafeTestsets

@safetestset "HMatrix" begin
    using HMatrices
    using Clusters
    using HMatrices: ACA, PartialACA
    using ComputationalResources
    using LinearAlgebra
    N    = 1000
    data = rand(Clusters.Point{2,Float64},N)
    splitter   = Clusters.CardinalitySplitter(nmax=128)
    clt  = Clusters.ClusterTree(data,splitter;reorder=true)
    adm  = Clusters.WeakAdmissibilityStd()
    bclt = Clusters.BlockTree(clt,clt,adm)
    f(x,y)::ComplexF64 = x==y ? 0.0 : exp(im*LinearAlgebra.norm(x-y))/LinearAlgebra.norm(x-y)
    M    = LazyMatrix(f,data,data)
    comp = HMatrices.ACA(rtol=1e-6)
    @testset "Assembly CPU1" begin
        H    = HMatrix(M,bclt,comp) |> Matrix
        @test norm(H-M,2) < comp.rtol*norm(M)
        @test norm(H-M,2) < comp.rtol*norm(M)
    end
    @testset "Assembly CPUThreads" begin
        H    = HMatrix(CPUThreads(),M,bclt,comp) |> Matrix
        @test norm(H-M,2) < comp.rtol*norm(M)
        @test norm(H-M,2) < comp.rtol*norm(M)
    end
end
