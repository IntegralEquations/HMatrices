using SafeTestsets

@safetestset "LazyMatrix" begin
    using HierarchicalMatrices
    using HierarchicalMatrices: LazyMatrix
    @testset "Simple ops" begin
        let
            m,n    = 10,5
            f(x,y) = x+y
            X      = rand(m)
            Y      = rand(n)
            M      = [f(x,y) for x in X, y in Y]
            lz     = LazyMatrix(f,X,Y)
            @test lz == M
            @test lz[7,3] == M[7,3]
            @test lz[:,2] == M[:,2]
            @test lz[end,:] == M[end,:]
        end
    end
end
