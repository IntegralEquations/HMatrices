using SafeTestsets

@safetestset "FlexMatrix" begin
    using HierarchicalMatrices.Matrices
    using HierarchicalMatrices.Matrices: num_elements, compression_rate, pushcol!
    @testset "Simple ops" begin
        let
            m,n    = 10,5
            data   = [rand(m) for _=1:n]
            F      = FlexMatrix(data)
            @test size(F) == (m,n)
            @test all(F[:,k] ==  data[k] for k=1:n)
            M = Matrix(F)
            @test F == FlexMatrix(M)
            pushcol!(F,rand(m),rand(m))
            @test size(F) == (10,n+2)
            push!(data,rand(m+1)) #invalidate data as a "matrix"
            @test_throws AssertionError FlexMatrix(data)
            @test_skip pushcol!(F,rand(m+1)) #this should throw an exception
        end
    end
end
