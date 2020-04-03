################################################################################
## LowerTriangular
################################################################################

const HUnitLowerTriangular = UnitLowerTriangular{<:Any,<:AbstractHMatrix}
const HLowerTriangular     = LowerTriangular{<:Any,<:AbstractHMatrix}
const HUpperTriangular     = UpperTriangular{<:Any,<:AbstractHMatrix}
const HUnitUpperTriangular = UnitUpperTriangular{<:Any,<:AbstractHMatrix}

getdata(L::HUnitLowerTriangular) = hasdata(L) ? UnitLowerTriangular(getdata(L.data)) : ()
hasdata(L::HUnitLowerTriangular) = hasdata(L.data)
colrange(L::HUnitLowerTriangular) = colrange(L.data)
rowrange(L::HUnitLowerTriangular) = rowrange(L.data)
isleaf(L::HUnitLowerTriangular)  = isleaf(L.data)
pivot(L::HUnitLowerTriangular)  = pivot(L.data)
blocksize(L::HUnitLowerTriangular)  = blocksize(L.data)
function getblock(L::HUnitLowerTriangular,i,j)
    chd = getblock(L.data,i,j)
    i===j ? UnitLowerTriangular(chd) : i>j ? chd : zero(chd)
end
Base.getindex(L::HUnitLowerTriangular,block::BlockIndex) = getblock(L,block.indices...)

getdata(L::HLowerTriangular) = hasdata(L) ? LowerTriangular(getdata(L.data)) : ()
hasdata(L::HLowerTriangular) = hasdata(L.data)
colrange(L::HLowerTriangular) = colrange(L.data)
rowrange(L::HLowerTriangular) = rowrange(L.data)
isleaf(L::HLowerTriangular)  = isleaf(L.data)
pivot(L::HLowerTriangular)  = pivot(L.data)
blocksize(L::HLowerTriangular)  = blocksize(L.data)
function getblock(L::HLowerTriangular,i,j)
    chd = getblock(L.data,i,j)
    i===j ? LowerTriangular(chd) : i>j ? chd : zero(chd)
end
Base.getindex(L::HLowerTriangular,block::BlockIndex) = getblock(L,block.indices...)

function ldiv!(L::Union{HUnitLowerTriangular,HLowerTriangular},b::AbstractMatrix)
    if isleaf(L)
        data = getdata(L)
        ldiv!(data,b) # b <-- L\b
    else
        @assert !hasdata(L) #only leaves are allowed to have data for the inversion
        shift = pivot(L) .- 1
        m,n = blocksize(L)
        @assert m === n
        for i=1:m
            irows  = colrange(L[block(i,i)]) .- shift[2]
            bi     = view(b,irows,:)
            for j=1:(i-1)# j<i
                jrows  = colrange(L[block(i,j)]) .- shift[2]
                bj     = view(b,jrows,:)
                mul!(bi,L[block(i,j)],bj,-1,true)
            end
            ldiv!(L[block(i,i)],bi)
        end
    end
    return b
end

ldiv!(L::Union{HUnitLowerTriangular,HLowerTriangular},R::RkMatrix) = (ldiv!(L,R.A); R)

function ldiv!(L::Union{HUnitLowerTriangular,HLowerTriangular},X::AbstractHMatrix)
    if isleaf(X)
        data = getdata(X)
        ldiv!(L,data)
    elseif isleaf(L) # X not a leaf, but L is a leaf
        @error "you probably should not be here..."
        # setdata!(H,Matrix(H))
        # H.children = ()
        # data = getdata(H)
        # ldiv!(L,data)
    else
        m,n = blocksize(L)
        @assert m === n
        for k=1:blocksize(X,2)
            for i=1:m
                for j=1:(i-1)# j<i
                    mul!(X[block(i,k)],L[block(i,j)],X[block(j,k)],-1,true)
                    flush_tree!(X[block(i,k)])
                end
                ldiv!(L[block(i,i)],X[block(i,k)])
            end
        end
    end
    return X
end

function rdiv!(b::StridedMatrix,L::Union{HUnitLowerTriangular,HLowerTriangular})
    if isleaf(L)
        data = getdata(L)
        rdiv!(b,data) # b <-- b/L
    else
        @assert !hasdata(L) #only leaves are allowed to have data for the inversion
        shift = pivot(L) .- 1 |> reverse
        m,n   = blocksize(L)
        @assert m === n
        for i=m:-1:1
            icols  = rowrange(L[block(i,i)]) .- shift[1]
            bi     = view(b,:,icols)
            for j=m:-1:(i+1)
                jcols  = rowrange(L[block(j,i)]) .- shift[1]
                bj     = view(b,:,jcols)
                mul!(bi,bj,L[block(j,i)],-1,true)
            end
            rdiv!(bi,L[block(i,i)])
        end
    end
    return b
end

# FIXME: rewrite the function below to avoid the allocations
function rdiv!(R::RkMatrix,L::Union{HUnitLowerTriangular,HLowerTriangular})
    Bt = rdiv!(Matrix(R.Bt),L)
    # R.B = adjoint(Bt) |> Matrix
    adjoint!(R.B,Bt)
    return R
end

function rdiv!(X::AbstractHMatrix,L::Union{HUnitLowerTriangular,HLowerTriangular})
    if isleaf(X)
        data = getdata(X)
        rdiv!(data,L) # b <-- b/L
    elseif isleaf(L)
        @error "should not happen"
    else
        @assert !hasdata(X) #only leaves are allowed to have data for the inversion
        m,n   = blocksize(L)
        @assert m === n
        for k=1:blocksize(X,1)
            for i=m:-1:1
                for j=m:-1:(i+1)
                    mul!(X[block(k,i)],X[block(k,j)],L[block(j,i)],-1,true)
                    flush_tree!(X[block(k,i)])
                end
                rdiv!(X[block(k,i)],L[block(i,i)])
            end
        end
    end
    return X
end

################################################################################
## UpperTriangular
################################################################################

getdata(L::HUnitUpperTriangular) = hasdata(L) ? UnitUpperTriangular(getdata(L.data)) : ()
hasdata(L::HUnitUpperTriangular) = hasdata(L.data)
colrange(L::HUnitUpperTriangular) = colrange(L.data)
rowrange(L::HUnitUpperTriangular) = rowrange(L.data)
isleaf(L::HUnitUpperTriangular)  = isleaf(L.data)
pivot(L::HUnitUpperTriangular)  = pivot(L.data)
blocksize(L::HUnitUpperTriangular)  = blocksize(L.data)
function getblock(L::HUnitUpperTriangular,i,j)
    chd = getblock(L.data,i,j)
    i===j ? UnitUpperTriangular(chd) : i<j ? chd : zero(chd)
end
Base.getindex(L::HUnitUpperTriangular,block::BlockIndex) = getblock(L,block.indices...)

getdata(L::HUpperTriangular) = hasdata(L) ? UpperTriangular(getdata(L.data)) : ()
hasdata(L::HUpperTriangular) = hasdata(L.data)
colrange(L::HUpperTriangular) = colrange(L.data)
rowrange(L::HUpperTriangular) = rowrange(L.data)
isleaf(L::HUpperTriangular)  = isleaf(L.data)
pivot(L::HUpperTriangular)  = pivot(L.data)
blocksize(L::HUpperTriangular)  = blocksize(L.data)
function getblock(L::HUpperTriangular,i,j)
    chd = getblock(L.data,i,j)
    i===j ? UpperTriangular(chd) : i<j ? chd : zero(chd)
end
Base.getindex(L::HUpperTriangular,block::BlockIndex) = getblock(L,block.indices...)

function ldiv!(U::Union{HUpperTriangular,HUnitUpperTriangular},b::AbstractMatrix)
    if isleaf(U)
        data = getdata(U)
        ldiv!(data,b) # b <-- U\b
    else
        @assert !hasdata(U) #only leaves are allowed to have data for the inversion
        shift = pivot(U) .- 1
        m,n = blocksize(U)
        @assert m === n
        for i=m:-1:1
            irows  = colrange(U[block(i,i)]) .- shift[2]
            bi     = view(b,irows,:)
            for j=i+1:n # j>i
                jrows  = colrange(U[block(i,j)]) .- shift[2]
                bj     = view(b,jrows,:)
                mul!(bi,U[block(i,j)],bj,-1,true)
            end
            ldiv!(U[block(i,i)],bi)
        end
    end
    return b
end

ldiv!(U::Union{HUpperTriangular,HUnitUpperTriangular},R::RkMatrix)     = (ldiv!(U,R.A); R)

function ldiv!(U::Union{HUpperTriangular,HUnitUpperTriangular},X::AbstractHMatrix)
    if isleaf(X)
        data = getdata(X)
        ldiv!(U,data)
    elseif isleaf(U) # X not a leaf, but U is a leaf
        @error "you should not have entered here..."
        # setdata!(H,Matrix(H))
        # H.children = ()
        # data = getdata(H)
        # ldiv!(U,data)
    else
        m,n = blocksize(U)
        @assert m === n
        for k=1:blocksize(X,2)
            for i in reverse(1:m)
                for j=i+1:m #j>i
                    mul!(X[block(i,k)],U[block(i,j)],X[block(j,k)],-1,true)
                    flush_tree!(X[block(i,k)])
                end
                ldiv!(U[block(i,i)],X[block(i,k)])
            end
        end
    end
    return X
end

function rdiv!(b::StridedMatrix,U::Union{HUnitUpperTriangular,HUpperTriangular})
    if isleaf(U)
        data = getdata(U)
        rdiv!(b,data) # b <-- b/L
    else
        @assert !hasdata(U) #only leaves are allowed to have data for the inversion
        shift = pivot(U) .- 1 |> reverse
        m,n   = blocksize(U)
        @assert m === n
        for i=1:m
            icols  = rowrange(U[block(i,i)]) .- shift[1]
            bi     = view(b,:,icols)
            for j=1:(i-1)
                jcols  = rowrange(U[block(j,i)]) .- shift[1]
                bj     = view(b,:,jcols)
                mul!(bi,bj,U[block(j,i)],-1,true)
            end
            rdiv!(bi,U[block(i,i)])
        end
    end
    return b
end

function rdiv!(R::RkMatrix,U::Union{HUnitUpperTriangular,HUpperTriangular})
    Bt = rdiv!(Matrix(R.Bt),U)
    adjoint!(R.B,Bt)
    return R
end

function rdiv!(X::AbstractHMatrix,U::Union{HUnitUpperTriangular,HUpperTriangular})
    if isleaf(X)
        data = getdata(X)
        rdiv!(data,U) # b <-- b/L
    elseif isleaf(U)
        @error "should not happen"
    else
        @assert !hasdata(U) #only leaves are allowed to have data for the inversion
        m,n   = blocksize(U)
        @assert m === n
        for k=1:blocksize(X,1)
            for i=1:m
                for j=1:(i-1)
                    mul!(X[block(k,i)],X[block(k,j)],U[block(j,i)],-1,true)
                end
                rdiv!(X[block(k,i)],U[block(i,i)])
            end
        end
    end
    return X
end
