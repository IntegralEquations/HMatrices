"""
    (+)(A::Union{Matrix,RkMatrix,HMatrix},B::Union{Matrix,RkMatrix,HMatrix}) --> C

Two argument addition. When operating on `Union{Matrix,RkMatrix,HMatrix}`, the result `C` is returned in the *natural format*, as
described in the table below:

| `(+)(A,B)` | `B::Matrix`  | `B::RkMatrix`   | `B::HMatrix` |
|:-----|:---:|:---:|:---:|
|`A::Matrix`  | `C::Matrix`  | `C::Matrix` | `C::HMatrix` |
|`A::RkMatrix`  | `C::Matrix` | `C::RkMatrix` | `C::HMatrix` |
|`A::HMatrix`  | `C::HMatrix` | `C::HMatrix` | `C::HMatrix` |
"""

#1.2
(+)(M::Matrix,R::RkMatrix) = M + Matrix(R)

#1.3
(+)(M::Matrix,H::HMatrix)  = axpby!(true,M,true,deepcopy(H))
(-)(M::Matrix,H::HMatrix)  = axpby!(-,M,true,deepcopy(H))

#2.1
(+)(R::RkMatrix,M::Matrix)  = (+)(M,R)
(-)(R::RkMatrix,M::Matrix)  = (-)(M,R)

#2.2
function (+)(R::RkMatrix,S::RkMatrix)
    Anew  = hcat(R.A,S.A)
    Bnew  = hcat(R.B,S.B)
    return RkMatrix(Anew,Bnew)
end
(-)(R::RkMatrix,S::RkMatrix) = (+)(R::RkMatrix,rmul!(S,-1))

#2.3
(+)(R::RkMatrix,H::HMatrix) = axpby!(true,R,true,deepcopy(H))
(-)(R::RkMatrix,H::HMatrix) = axpby!(-1,R,true,deepcopy(H))

#3.1
(+)(H::HMatrix,M::Matrix) = (+)(M,H)
(-)(H::HMatrix,M::Matrix) = (-)(M,H)

#3.2
(+)(H::HMatrix,R::RkMatrix) = (+)(R,H)
(-)(H::HMatrix,R::RkMatrix) = (-)(R,H)

#3.3
(+)(H::HMatrix,S::HMatrix) = axpby!(true,H,true,deepcopy(S))
(-)(H::HMatrix,S::HMatrix) = axpby!(-1,H,true,deepcopy(S))

# Below we define the more general method axpby!, and then base the other addition operations on it

#1.2
axpby!(a,X::Matrix,b,Y::RkMatrix,compress=identity) = axpby!(a,RkMatrix(X),b,Y,compress)

#1.3
function axpby!(a,X::Matrix,b,Y::HMatrix,compress=identity)
    rmul!(Y,b)
    if hasdata(Y)
        if getdata(Y) isa Matrix
            axpby!(a,X,true,getdata(Y))
        else
            axpby!(a,X,true,getdata(Y),compress)
        end
    else
        isleaf(Y) || warn("adding data to an internal node. This may not be what you intended. ")
        data = a*X
        setdata!(Y,data)
    end
    return Y
end

#2.1
axpby!(a,X::RkMatrix,b,Y::Matrix,compress=identity) = axpby!(a,Matrix(X),b,Y,compress=identity)

#2.2
function axpby!(a,X::RkMatrix,b,Y::RkMatrix,compress=identity)
    rmul!(Y,b)
    m,n = size(X)
    if m<n
        Y.A   = hcat(a*X.A,Y.A) |> compress
        Y.B   = hcat(X.B,Y.B)   |> compress
    else
        Y.A   = hcat(X.A,Y.A)   |> compress
        Y.B   = hcat(a*X.B,Y.B) |> compress
    end
    return Y
end

#2.3
function axpby!(a,X::RkMatrix,b,Y::HMatrix,compress=identity)
    rmul!(Y,b)
    if hasdata(Y)
        data = getdata(Y)
        if data isa Matrix
            axpby!(a,X,true,getdata(Y))
        else
            axpby!(a,X,true,getdata(Y),compress)
        end
    else
        data = a*X
        setdata!(Y,data)
    end
    return Y
end

#3.1
function axpby!(a,X::HMatrix,b,Y::Matrix)
    rmul!(Y,b)
    shift = pivot(X) .- 1
    for block in PreOrderDFS(X)
        irange = rowrange(block) .- shift[1]
        jrange = colrange(block) .- shift[2]
        if hasdata(block)
            axpby!(a,getdata(block),true,view(Y,irange,jrange))
        end
    end
    return Y
end

#3.2
function axpby!(a,X::HMatrix,b,Y::RkMatrix,compress=identity)
    R = RkMatrix(X)
    axpby!(a,R,b,Y,compress)
end

function axpby!(a,X::HMatrix,b,Y::HMatrix,compress)
    rmul!(Y,b)
    if hasdata(X)
        axpby!(a,getdata(X),true,Y,compress)
    end
    for (bx,by) in zip(getchildren(X),getchildren(Y))
        axpby!(a,bx,true,by,compress)
    end
    return Y
end

# some extra cases
function axpby!(a,X::UniformScaling,b,Y::HMatrix)
    rmul!(Y,b)
    if hasdata(Y)
        data = getdata(Y)
        @assert data isa Matrix
        n = min(size(data)...)
        for i=1:n
            data[i,i] += a*X.λ
        end
    else
        n = min(blocksize(Y)...)
        for i=1:n
            axpby!(a,X,true,getchildren(Y,i,i))
        end
    end
    return Y
end

@inline function axpby!(a,X::AbstractSparseArray{<:Any,<:Any,2},b,Y::HMatrix)
    rmul!(Y,b)
    if isleaf(Y)
        Xblock = X[rowrange(Y),colrange(Y)] #extract sparse block
        if !iszero(Xblock)
            data = getdata(Y)
            axpy!(a,Xblock,data)
        end
        return Y
    else
        for child in getchildren(Y)
            axpby!(a,X,true,child)
        end
    end
    return Y
end
axpy!(a,X::AbstractSparseArray,Y::HMatrix) = axpby!(a,X::AbstractSparseArray,true,Y::HMatrix)

(+)(X::AbstractSparseArray{<:Any,<:Any,2},Y::HMatrix) = axpy!(true,X,deepcopy(Y))
(+)(X::HMatrix,Y::AbstractSparseArray{<:Any,<:Any,2}) = Y+X

function axpy!(a,X::AbstractSparseArray{<:Any,<:Any,2},Y::RkMatrix)
    M = Matrix(Y)
    axpy!(a,X,M)
    R = RkMatrix(M)
    Y.A = R.A
    Y.B = R.B
    return Y
end

(+)(X::UniformScaling,Y::HMatrix) = axpby!(true,X,true,deepcopy(Y))
(+)(X::HMatrix,Y::UniformScaling) = Y+X
