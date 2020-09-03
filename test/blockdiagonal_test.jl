using Test
using SafeTestsets

@safetestset "BlockDiagonal" begin
    using HMatrices
    using HMatrices: BlockDiagonal
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
    f(x,y)::ComplexF64 = x==y ? sum(x.-y) : exp(im*LinearAlgebra.norm(x.-y))/LinearAlgebra.norm(x.-y)
    M    = LazyMatrix(f,data,data)
    comp = HMatrices.PartialACA(rtol=1e-6)
    blocks = [rand(10,10) for  _ in 1:10]
    B    = BlockDiagonal(blocks)
    @testset "Basics" begin
        @test B[1:10,1:10]   == blocks[1]
        @test B[11:20,11:20] == blocks[2]
        Binv = inv(B)
        @test Binv[1:10,1:10] == inv(blocks[1])
        @test Binv[11:20,11:20]  == inv(blocks[2])
    end
    @testset "Multiplication" begin
        x    = rand(size(B,2))
        @test B*x ≈ Matrix(B)*x
        y    = similar(x)
        @test mul!(y,B,x) ≈ B*x
    end
end
