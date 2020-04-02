using SafeTestsets

@safetestset "Multiplication" begin
    using HierarchicalMatrices
    using HierarchicalMatrices: ACA, PartialACA, RkMatrix
    using LinearAlgebra

    N,r    = 500, 3
    T      = ComplexF64
    data   = HierarchicalMatrices.points_on_cylinder(N,1,3/sqrt(N))
    splitter   = Clusters.GeometricMinimalSplitter()
    clt  = Clusters.ClusterTree(data,splitter;reorder=true)
    adm  = Clusters.StrongAdmissibilityStd(100)
    bclt = Clusters.BlockTree(clt,clt,adm)
    comp = HierarchicalMatrices.PartialACA(rtol=1e-6)
    f(x,y)::ComplexF64 = x==y ? 0.0 : exp(im*LinearAlgebra.norm(x-y))/LinearAlgebra.norm(x-y)
    L    = LazyMatrix(f,data,data)
    H    = HMatrix(L,bclt,comp)
    _H   = Matrix(H)
    R    = RkMatrix(rand(T,N,r),rand(T,N,r))
    _R   = Matrix(R)
    M    = rand(T,N,N)
    _M   = Matrix(M)
    a    = rand(); b = rand();

    @testset "gemv" begin
        @testset "CPU1 mul!" begin
            x    = rand(ComplexF64,N)
            y    = similar(x)
            @test mul!(y,R,x) ≈ _R*x
            @test mul!(y,adjoint(R),x) ≈ adjoint(_R)*x
            @test mul!(y,H,x) ≈ _H*x
            @test mul!(y,adjoint(H),x) ≈ adjoint(_H)*x
            tmp = b*y + a*_H*x
            @test mul!(y,H,x,a,b) ≈ tmp
        end
        @testset "CPUThreads mul!" begin
            x    = rand(ComplexF64,N)
            y    = zero(x)
            @test mul!(CPUThreads(),y,H,x,true,false) ≈ _H*x
            tmp = b*y + a*_H*x
            @test mul!(CPUThreads(),y,H,x,a,b) ≈ tmp
        end
    end
    #

    @testset "gemm" begin
        @testset "*" begin
            @testset "1.2" begin
                tmp = M*R
                @test tmp isa RkMatrix
                @test M*R ≈ _M*_R
            end
            @testset "1.3" begin
                tmp = M*H
                @test tmp isa Matrix
                @test tmp ≈ _M*_H
            end
            @testset "2.1" begin
                tmp = R*M
                @test tmp isa RkMatrix
                @test tmp ≈ _R*_M
            end
            @testset "2.2" begin
                tmp = R*R
                @test tmp isa RkMatrix
                @test tmp ≈ _R*_R
            end
            @testset "2.3" begin
                tmp = R*H
                @test tmp isa RkMatrix
                @test tmp ≈ _R*_H
            end
            @testset "3.1" begin
                tmp = H*M
                @test tmp isa Matrix
                @test tmp ≈ _H*_M
            end
            @testset "3.2" begin
                tmp = H*R
                @test tmp isa RkMatrix
                @test tmp ≈ _H*_R
            end
            @testset "3.3" begin
                # tmp = H*H
                # @test tmp isa HMatrix
                # @test tmp ≈ _H*_H
            end
        end
        
        @testset "mul!" begin
            #cases whose target is a full matrix
            @testset "1.1.2" begin
                exact = b*_M + a*_M*_R
                tmp   = copy(M)
                @test mul!(tmp,M,R,a,b) ≈ exact
            end
            @testset "1.1.3" begin
                exact = b*_M + a*_M*_H
                tmp   = copy(M)
                @test mul!(tmp,M,H,a,b) ≈ exact
            end
            @testset "1.2.1" begin
                exact = b*_M + a*_R*_M
                tmp   = copy(M)
                @test mul!(tmp,R,M,a,b) ≈ exact
                exact = b*_M + a*adjoint(_R)*_M
                tmp   = copy(M)
                @test mul!(tmp,adjoint(R),M,a,b) ≈ exact
            end
            @testset "1.2.2" begin
                exact = b*_M + a*_R*_R
                tmp   = copy(M)
                @test mul!(tmp,R,R,a,b) ≈ exact
            end
            @testset "1.2.3" begin
                exact = b*_M + a*_R*_H
                tmp   = copy(M)
                @test mul!(tmp,R,H,a,b) ≈ exact
            end
            @testset "1.3.1" begin
                exact = b*_M + a*_H*_M
                tmp   = copy(M)
                @test mul!(tmp,H,M,a,b) ≈ exact
                exact = b*_M + a*adjoint(_H)*_M
                tmp   = copy(M)
                @test mul!(tmp,adjoint(H),M,a,b) ≈ exact
            end
            @testset "1.3.2" begin
                exact = b*_M + a*_H*_R
                tmp   = copy(M)
                @test mul!(tmp,H,R,a,b) ≈ exact
            end
            @testset "1.3.3" begin
                exact = b*_M + a*_H*_H
                tmp   = copy(M)
                @test mul!(tmp,H,H,a,b) ≈ exact
            end
            #cases whose target is a sparse matrix
            @testset "2.1.1" begin
                exact = b*_R + a*_M*_M
                tmp   = copy(R)
                @test mul!(tmp,M,M,a,b) ≈ exact
            end
            @testset "2.1.2" begin
                exact = b*_R + a*_M*_R
                tmp   = deepcopy(R)
                @test mul!(tmp,M,R,a,b) ≈ exact
            end
            @testset "2.1.3" begin
                exact = b*_R + a*_M*_H
                tmp   = deepcopy(R)
                @test mul!(tmp,M,H,a,b) ≈ exact
            end
            @testset "2.2.1" begin
                exact = b*_R + a*_R*_M
                tmp   = deepcopy(R)
                @test mul!(tmp,R,M,a,b) ≈ exact
            end
            @testset "2.2.2" begin
                exact = b*_R + a*_R*_R
                tmp   = deepcopy(R)
                @test mul!(tmp,R,R,a,b) ≈ exact
            end
            @testset "2.2.3" begin
                exact = b*_R + a*_R*_H
                tmp   = deepcopy(R)
                @test mul!(tmp,R,H,a,b) ≈ exact
            end
            @testset "2.3.1" begin
                exact = b*_R + a*_H*_M
                tmp   = deepcopy(R)
                @test mul!(tmp,H,M,a,b) ≈ exact
            end
            @testset "2.3.2" begin
                exact = b*_R + a*_H*_R
                tmp   = deepcopy(R)
                @test mul!(tmp,H,R,a,b) ≈ exact
            end
            @testset "2.3.2" begin
                exact = b*_R + a*_H*_R
                tmp   = deepcopy(R)
                @test mul!(tmp,H,R,a,b) ≈ exact
            end
            @testset "2.3.3" begin
                exact = b*_R + a*_H*_H
                tmp   = deepcopy(R)
                @test mul!(tmp,H,H,a,b) ≈ exact
            end
            #cases where target is a hierarchical matrix
            @testset "3.1.1" begin
                exact = b*_H + a*_M*_M
                tmp   = deepcopy(H)
                @test mul!(tmp,M,M,a,b) ≈ exact
            end
            @testset "3.1.2" begin
                exact = b*_H + a*_M*_R
                tmp   = deepcopy(H)
                @test mul!(tmp,M,R,a,b) ≈ exact
            end
            @testset "3.1.3" begin
                exact = b*_H + a*_M*_H
                tmp   = deepcopy(H)
                @test mul!(tmp,M,H,a,b) ≈ exact
            end
            @testset "3.2.1" begin
                exact = b*_H + a*_R*_M
                tmp   = deepcopy(H)
                @test mul!(tmp,R,M,a,b) ≈ exact
            end
            @testset "3.2.2" begin
                exact = b*_H + a*_R*_R
                tmp   = deepcopy(H)
                @test mul!(tmp,R,R,a,b) ≈ exact
            end
            @testset "3.2.3" begin
                exact = b*_H + a*_R*_H
                tmp   = deepcopy(H)
                @test mul!(tmp,R,H,a,b) ≈ exact
            end
            @testset "3.3.1" begin
                exact = b*_H + a*_H*_M
                tmp   = deepcopy(H)
                @test mul!(tmp,H,M,a,b) ≈ exact
            end
            @testset "3.3.2" begin
                exact = b*_H + a*_H*_R
                tmp   = deepcopy(H)
                @test mul!(tmp,H,R,a,b) ≈ exact
            end
            @testset "3.3.3" begin
                exact = b*_H + a*_H*_H
                tmp   = deepcopy(H)
                @test mul!(tmp,H,H,a,b) ≈ exact
            end
        end#mul!
    end#gemm
    @testset "flush data" begin
        using HierarchicalMatrices: flush_to_children!, flush_to_leaves!, hasdata, setdata!, flush_tree!

        @testset "flush_to_children" begin
            exact = _H + _R
            tmp   = deepcopy(H)
            setdata!(tmp,R)
            tmp ≈ exact
            @test hasdata(tmp) == true
            HierarchicalMatrices.flush_to_children!(tmp)
            @test hasdata(tmp) == false
            @test tmp ≈ exact
        end
        @testset "flush_to_leaves" begin
            exact = _H + _R
            tmp   = deepcopy(H)
            setdata!(tmp,R)
            tmp ≈ exact
            @test hasdata(tmp) == true
            HierarchicalMatrices.flush_to_leaves!(tmp)
            @test hasdata(tmp) == false
            @test tmp ≈ exact
        end
        @testset "flush_tree" begin
            exact = _H + _R
            tmp   = deepcopy(H)
            setdata!(tmp,R)
            tmp ≈ exact
            @test hasdata(tmp) == true
            HierarchicalMatrices.flush_tree!(tmp)
            @test hasdata(tmp) == false
            @test tmp ≈ exact
        end
    end
end
