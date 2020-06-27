LU(H::HMatrix) = LU(H,Int[],0)

const HLU = LU{<:Any,<:HMatrix}

unit_tril(H::HLU) = unit_tril!(deepcopy(H.factors))
triu(H::HLU)  = triu!(deepcopy(H.factors))

function unit_tril!(H::HMatrix)
    T = eltype(H)
    m,n = blocksize(H)
    if hasdata(H)
        data = (getdata(H)) |> UnitLowerTriangular |> full!
        setdata!(H,data)
    end
    for i=1:m
        for j=1:n
            b = H[block(i,j)]
            if j>i
                setdata!(b,zero(RkMatrix{T},size(b)...))
                setchildren!(b,())
                b.admissible=true
            elseif j<i
                # do nothing
            else
                unit_tril!(b)
            end
        end
    end
    return H
end

function triu!(H::HMatrix)
    T = eltype(H)
    m,n = blocksize(H)
    if hasdata(H)
        data = (getdata(H)) |> UpperTriangular |> full!
        setdata!(H,data)
    end
    for i=1:m
        for j=1:n
            b = H[block(i,j)]
            if j<i
                setdata!(b,zero(RkMatrix{T},size(b)...))
                setchildren!(b,())
                b.admissible=true
            elseif j>i
                # do nothing
            else
                triu!(b)
            end
        end
    end
    return H
end

function Base.getproperty(LU::HLU,s::Symbol)
    T = eltype(LU)
    if s == :L
        return unit_tril(LU)
    elseif s == :U
        return triu(LU)
    else
        return getfield(LU,s)
    end
end

function lu!(M::HMatrix)
    #perform the lu decomposition of M in place
    _lu!(M)
    #wrap the result in the LU structure
    return LU(M)
end

function _lu!(M::HMatrix)
    if isleaf(M)
        data = getdata(M)
        @assert data isa Matrix
        lu!(data,Val(false))#Val(false) for pivot of dense lu factorization
    else
        @assert !hasdata(M)
        m,n = blocksize(M)
        for i=1:m
            _lu!(M[block(i,i)])
            for j=i+1:n
                ldiv!(UnitLowerTriangular(M[block(i,i)]),M[block(i,j)])
                rdiv!(M[block(j,i)],UpperTriangular(M[block(i,i)]))
                mul!(M[block(j,j)],M[block(j,i)],M[block(i,j)],-1,1)
                flush_tree!(M[block(j,j)])
            end
        end
    end
    return M
end

lu(M::HMatrix) = lu!(deepcopy(M))
