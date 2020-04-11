using SafeTestsets

@safetestset "RkMatrix" begin
    using HierarchicalMatrices
    using HierarchicalMatrices: FlexMatrix, RkFlexMatrix, RkMatrix, num_elements, compression_rate
    @testset "Simple ops" begin
        let
            m,n,r = 10,20,5
            A = rand(m,r)
            B = rand(n,r)
            M = A*adjoint(B)
            Rk = RkMatrix(A,B)
            @test size(rand(RkMatrix{ComplexF64},m,n,r)) == (m,n)
            @test eltype(rand(RkMatrix,m,n,r)) == Float64
            @test Rk ≈ M
            @test Rk.Bt == adjoint(B)
            @test Rk.At == adjoint(A)
            @test Rk[1:3,4:10] ≈ M[1:3,4:10]
            @test hcat(Rk,Rk) ≈ hcat(M,M)
            @test vcat(Rk,Rk) ≈ vcat(M,M)
            @test num_elements(Rk) == (m+n)*r
            @test compression_rate(Rk) == (m+n)*r / (m*n)
        end
    end
end
