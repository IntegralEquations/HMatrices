const HUnitLowerTriangular       = UnitLowerTriangular{<:Any,<:AbstractHMatrix}
const HLowerTriangular           = LowerTriangular{<:Any,<:AbstractHMatrix}

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
                    flush_tree!(X[block[i,k]])
                end
                ldiv!(UnitLowerTriangular(L[block(i,i)]),X[block(i,k)])
            end
        end
    end
    return X
end

################################################################################
## UpperTriangular
################################################################################
const HUpperTriangular       = UpperTriangular{<:Any,<:AbstractHMatrix}

getdata(L::UpperTriangular) = hasdata(L) ? UpperTriangular(getdata(L.data)) : ()
hasdata(L::UpperTriangular) = hasdata(L.data)
colrange(L::UpperTriangular) = colrange(L.data)
rowrange(L::UpperTriangular) = rowrange(L.data)
isleaf(L::UpperTriangular)  = isleaf(L.data)
pivot(L::UpperTriangular)  = pivot(L.data)
blocksize(L::UpperTriangular)  = blocksize(L.data)
function getblock(L::UpperTriangular,i,j)
    chd = getblock(L.data,i,j)
    i===j ? UpperTriangular(chd) : i<j ? chd : zero(chd)
end
Base.getindex(L::UpperTriangular,block::BlockIndex) = getblock(L,block.indices...)

function ldiv!(U::HUpperTriangular,b::AbstractMatrix)
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

ldiv!(U::HUpperTriangular,R::RkMatrix)     = (ldiv!(U,R.A); R)

function ldiv!(U::HUpperTriangular,X::AbstractHMatrix)
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
                    flush_tree!(X[block[i,k]])
                end
                ldiv!(UpperTriangular(U[block(i,i)]),X[block(i,k)])
            end
        end
    end
    return X
end
