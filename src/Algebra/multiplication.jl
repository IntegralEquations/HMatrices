function LinearAlgebra.mul!(C::AbstractVector,H::HMatrix,F::AbstractVector,a::Number,b::Number)
    LinearAlgebra.rmul!(C,b)
    shift = _idx_pivot(H) .- 1
    for block in AbstractTrees.Leaves(H)
        irange = rowrange(block) .- shift[1]
        jrange = colrange(block) .- shift[2]
        data = getdata(block)
        LinearAlgebra.mul!(view(C,irange,:),data,view(F,jrange,:),a,1)
    end
    return C
end

function LinearAlgebra.mul!(C::AbstractVector,Rk::RkFlexMatrix,F::AbstractVector,a::Number,b::Number)
    LinearAlgebra.mul!(C,Rk.A,Rk.Bt*F,a,b)
end

function LinearAlgebra.mul!(C::Vector,Rk::RkMatrix,F::Vector,a,b)
    @info "here"
    LinearAlgebra.mul!(C,Rk.A,Rk.Bt*F,a,b)
end
