using SafeTestsets

@safetestset "RkFlexMatrix" begin
    using HierarchicalMatrices
    using HierarchicalMatrices: FlexMatrix, RkFlexMatrix, num_elements, compression_rate
    @testset "Simple ops" begin
        let
            T = ComplexF64
            m,n,r = 10,20,5
            A = FlexMatrix([rand(T,m) for _=1:r])
            B = FlexMatrix([rand(T,n) for _=1:r])
            M = A*adjoint(B)
            Rk = RkFlexMatrix(A,B)
            @test size(Rk) == (m,n)
            @test eltype(Rk) == T
            @test Rk ≈ M
            @test Rk.Bt == adjoint(B)
            @test Rk.At == adjoint(A)
            @test hcat(Rk,Rk) ≈ hcat(M,M)
            @test vcat(Rk,Rk) ≈ vcat(M,M)
        end
    end
end
