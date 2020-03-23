using SafeTestsets

@safetestset "Multiplication" begin
    using HierarchicalMatrices
    using HierarchicalMatrices: ACA, PartialACA, hmul!
    using LinearAlgebra
    @testset "CPU1 mul!" begin
        let
            N    = 1000
            data = rand(Geometry.Point{2,Float64},N)
            splitter   = Clusters.GeometricMinimalSplitter()
            clt  = Clusters.ClusterTree(data,splitter;reorder=true)
            bclt = Clusters.BlockTree(clt,clt)
            comp = HierarchicalMatrices.PartialACA(rtol=1e-4)
            f(x,y)::ComplexF64 = x==y ? 0.0 : exp(im*LinearAlgebra.norm(x-y))/LinearAlgebra.norm(x-y)
            M    = LazyMatrix(f,data,data)
            H    = HMatrix(M,bclt,comp)
            x    = rand(ComplexF64,N)
            y    = similar(x)
            # 3-arg mul
            mul!(y,H,x)
            @test norm(y-M*x) < 10*comp.rtol*norm(M)
            a = 1; b=2;
            tmp = b*y + a*M*x
            # 5-arg mul
            hmul!(y,H,x,a,b)
            @test norm(y-tmp) < 10*comp.rtol*norm(M)
        end
    end
    @testset "CPUThreads mul!" begin
        let
            N    = 1000
            data = Geometry.points_on_cylinder(N,1,3/sqrt(N))
            # data = rand(Geometry.Point{2,Float64},N)
            splitter   = Clusters.CardinalitySplitter()
            clt  = Clusters.ClusterTree(data,splitter;reorder=true)
            bclt = Clusters.BlockTree(clt,clt)
            comp = HierarchicalMatrices.PartialACA(rtol=1e-4)
            f(x,y)::ComplexF64 = x==y ? 0.0 : exp(im*LinearAlgebra.norm(x-y))/LinearAlgebra.norm(x-y)
            M    = LazyMatrix(f,data,data)
            H    = HMatrix(M,bclt,comp)
            x    = rand(ComplexF64,N)
            y    = zero(x)
            # 3-arg mul
            mul!(CPUThreads(),y,H,x)
            @test norm(y-M*x) < 10*comp.rtol*norm(M)
            # 5-arg mul
            a = 1; b=2;
            tmp = b*y + a*M*x
            hmul!(CPUThreads(),y,H,x,a,b)
            @test norm(y-tmp) < 10*comp.rtol*norm(M)
        end
    end
    @testset "plan_mul!" begin
        let
            # using HScheduler
            # N    = 4000
            # data = rand(Point{2,Float64},N)
            # splitter   = Clusters.CardinalitySplitter(nmax=128)
            # clt = ClusterTree(data,splitter;reorder=true)
            # bclt = BlockTree(clt,clt)
            # comp = HierarchicalMatrices.PartialACA(rtol=1e-6)
            # f(x,y)::ComplexF64 = x==y ? 0.0 : exp(im*LinearAlgebra.norm(x-y))/LinearAlgebra.norm(x-y)
            # M    = LazyMatrix(f,data,data)
            # H    = HMatrix(M,bclt,comp)
            # x    = rand(ComplexF64,N)
            # y  = zero(x)
            # tg = plan_mul!(y,H,x,1,0)
            # t1 = @elapsed HScheduler.execute.(tg.vertices)
            
            # # @test norm(y-M*x) < comp.rtol*norm(M)
            # tg = plan_mul!(y,H,x,1,0)
            # fill!(y,0)
            # # t2 = @elapsed HScheduler.run_parallel(tg)
            # t3 = @elapsed mul!(y,H,x,1,0,Val(false))
            # t4 = @elapsed mul!(y,H,x,1,0,Val(true))
            # @show t1, t3, t4
            # @test norm(y-M*x) < comp.rtol*norm(M)
        end
    end
end
