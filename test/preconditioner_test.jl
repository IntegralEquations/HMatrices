using Test
using SafeTestsets

@safetestset "Preconditioners" begin
    using HMatrices
    using Clusters
    using HMatrices: ACA, PartialACA
    using ComputationalResources
    using LinearAlgebra
    N    = 1000
    data = [(rand(),rand()) for _ in 1:N]
    splitter   = Clusters.CardinalitySplitter(nmax=128)
    clt  = Clusters.ClusterTree(data,splitter;reorder=true)
    adm  = Clusters.WeakAdmissibilityStd()
    bclt = Clusters.BlockTree(clt,clt,adm)
    f(x,y)::ComplexF64 = x==y ? sum(x.+y) : exp(im*LinearAlgebra.norm(x.-y))/LinearAlgebra.norm(x.-y)
    M    = LazyMatrix(f,data,data)
    comp = HMatrices.PartialACA(rtol=1e-6)
    H    = HMatrix(M,bclt,comp)
    @testset "Diagonal" begin
        using HMatrices: DiagonalPreconditioner
        P = DiagonalPreconditioner(H)
        @test P.matrix  == Diagonal(M)
        @test P.inverse == inv(Diagonal(M))
    end
    @testset "Block diagonal" begin
        using HMatrices: block_diag, BlockDiagonalPreconditioner
        P  = BlockDiagonalPreconditioner(H);
        x  = rand(ComplexF64,size(P,2))
    end
end
