_default_rtol(M) = sqrt(eps(eltype(M)))

function aca_full(M; atol::Real=0, rtol::Real = atol>0 ? 0 : _default_rtol(M),  norm::Function=LinearAlgebra.norm)
    @debug rtol, atol
    T = eltype(M)
    m,n = size(M)
    R = RkFlexMatrix{T}(undef,m,n,0)
    er = Inf
    enorm = norm(M) #exact norm
    while er > max(atol,rtol*enorm)
        @debug er
        (i,j) = argmax(abs.(M)).I
        δ       = M[i,j]
        if δ == 0
            return R
        else
            col,row = cross(M,i,j)
            LinearAlgebra.rdiv!(row,δ)
            @debug size(R)
            Matrices.pushcross!(R,col,row)
            M  = M - col*adjoint(row)
            er = norm(col)*norm(row) #estimate of error in norm
        end
    end
    return R
end

function cross(M,i,j)
    col = M[:,j]
    row = M[i,:]
    return col, row
end

function aca_partial(M; atol::Real=0, rtol::Real = atol>0 ? 0 : sqrt(eps(eltype(M))),  norm::Function=LinearAlgebra.norm)
    T = Base.eltype(M,Int,Int)
    A = Vector{T}[]
    B = Vector{T}[]
    I = [true for i = 1:size(M,1)]
    J = [true for i = 1:size(M,2)]
    i = 1
    er = Inf
    while er > max(atol,rtol*enorm)
        I[i] = false  # remove index i from allowed row
        b    = M[i,:]
        _update_row!(b,A,B,i)
        j    = nextcol(b,J)
        δ    = b[j]
        if δ == 0
            i = findfirst(x->x==true,J)
            isnothing(i) && error("aca did not converge")
        else # δ != 0
            rdiv!(b,δ) # b <-- b/δ
            J[j] = false
            a    = M[:,j] - R[:,j] # compute a column
            append!(R, a, b)
            er = norm(a)*norm(b) # approximate error
            i = nextrow(a,I)
        end
    end
    conj!.(B)
    return RkMatrix(A,B)
end

function nextcol(col,J)
    out = -1
    val = -Inf
    for n in 1:length(J)
        J[n] || continue
        tmp = abs(col[n])
        tmp < val && continue
        out = n
        val = tmp
    end
    return out
end
nextrow(row,I) = nextcol(row,I)

function _update_row!(b,A,B,i)
    for n = eachindex(A)
        axpy!(-A[n][i],B[n],b)
    end
    return b
end

function _update_col!(a,A,B,j)
    for (n,i) in enumerate(irange)
        a[n] = f(i,j+jshift) # evaluate col of original matrix
    end
    return a
end
