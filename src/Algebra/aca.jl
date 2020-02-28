function aca_full(M::Matrix{T}; atol::Real=0, rtol::Real = atol>0 ? 0 : sqrt(eps(T)),  norm::Function=LinearAlgebra.norm) where {T}
    @debug rtol, atol
    A = Vector{T}[]
    B = Vector{T}[]
    er = Inf
    enorm = norm(M) #exact norm
    while er > max(atol,rtol*enorm)
        @debug er
        (i,j) = argmax(abs.(M)).I
        δ     = M[i,j]
        if δ == 0
            return RkMatrix(A,B)
        else
            new_col = M[:,j]
            new_row = M[i,:]/δ
            push!(A,new_col)
            push!(B,new_row)
            M  = M - new_col*adjoint(new_row)
            er = norm(new_col)*norm(new_row) #estimate of error in norm
        end
    end
    return RkMatrix(A,B)
end
