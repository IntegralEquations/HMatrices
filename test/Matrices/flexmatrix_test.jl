using SafeTestsets

@safetestset "FlexMatrix" begin
    using HierarchicalMatrices.Matrices
    using HierarchicalMatrices.Matrices: num_elements, compression_rate
    @testset "Simple ops" begin
        let
            m,n    = 10,5
            data   = [rand(m) for _=1:n]
            F      = FlexMatrix(data)
            @test size(F) == (10,5)
            @test all(F[:,k] ==  data[k] for k=1:n)
            M = Matrix(F)
            @test F == FlexMatrix(M)
            push!(data,rand(m+1)) #invalidate data as a "matrix"
            @test_throws AssertionError FlexMatrix(data)
        end
    end
end
