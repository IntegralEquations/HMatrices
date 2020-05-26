################################################################################
## rmul!
################################################################################
function rmul!(R::RkMatrix,b::Number)
    m,n = size(R)
    if m>n
        rmul!(R.B,conj(b))
    else
        rmul!(R.A,b)
    end
    return R
end
(*)(a::Number,R::RkMatrix) = rmul!(deepcopy(R),a)
(*)(R::RkMatrix,a::Number) = rmul!(deepcopy(R),a)

@inline function rmul!(H::HMatrix,b::Number)
    b==true && (return H) # short circuit. If inlined, allows for compiler to optimize rmul!(H,true) to a no-op
    if hasdata(H)
        rmul!(getdata(H),b)
    end
    for child in getchildren(H)
        rmul!(child,b)
    end
    return H
end
(*)(a::Number,H::HMatrix) = rmul!(deepcopy(H),a)
(*)(H::HMatrix,a::Number) = rmul!(deepcopy(H),a)

################################################################################
## gemv
################################################################################

# R*x
function mul!(C::AbstractVector,Rk::RkMatrix,F::AbstractVector,a::Number,b::Number,buffer=similar(C,rank(Rk)))
    buffer = mul!(buffer,Rk.Bt,F)
    mul!(C,Rk.A,buffer,a,b)
    return C
end

# Rt*x
function mul!(y::AbstractVector,Rt::Adjoint{<:Any,<:RkMatrix},x::AbstractVector,a::Number,b::Number)
    R  = Rt.parent
    At = adjoint(R.A)
    buffer = At*x
    mul!(y,R.B,buffer,a,b)
    return y
end

# xt*R
const VecAdj{T} = Adjoint{T,Vector{T}}
function mul!(yt::VecAdj,xt::VecAdj,R::RkMatrix,a::Number,b::Number)
    mul!(yt.parent,adjoint(R),xt.parent,a,b)
    return yt
end

# CPU1
function mul!(resources::CPU1,C::AbstractVector,H::AbstractHMatrix,F::AbstractVector,a::Number=true,b::Number=false)
    rmul!(C,b)
    if hasdata(H)
        _multiply_leaf!(C,H,F,a,b)
    end
    for child in getchildren(H)
        mul!(C,child,F,a,true)
    end
    return C
end

function mul!(y::AbstractVector,adjH::Adjoint{<:Any,<:AbstractHMatrix},x::AbstractVector,a::Number=true,b::Number=false)
    rmul!(y,b)
    if isleaf(adjH)
        _multiply_leaf!(y,adjH,x,a,true)
    else
        for child in getchildren(adjH)
            mul!(y,child,x,a,true)
        end
    end
    return y
end

# CPUThreads
function mul!(resources::CPUThreads,C::AbstractVector,H::AbstractHMatrix,F::AbstractVector,a::Number=true,b::Number=false)
    rmul!(C,b)
    nthreads = Threads.nthreads()
    y        = [zero(C) for _ = 1:nthreads]
    @sync for leaf in Leaves(H)
        @spawn begin
            i = Threads.threadid()
            _multiply_leaf!(y[i],leaf,F,a,1)
        end
    end
    #reduction stage
    for i=1:nthreads
        axpy!(1,y[i],C)
    end
    return C
end

# default
mul!(C::AbstractVector,H::AbstractHMatrix,F::AbstractVector,a::Number,b::Number) = mul!(CPUThreads(),C,H,F,a,b)

function _multiply_leaf!(C::AbstractVector,H,F::AbstractVector,a,b)
    irange = rowrange(H)
    jrange = colrange(H)
    Cview  = view(C,irange)
    Fview  = view(F,jrange)
    data   = getdata(H)
    mul!(Cview,data,Fview,a,b)
    return C
end


################################################################################
## GEMM
################################################################################

"""
    (*)(A::Union{Matrix,RkMatrix,HMatrix},B::Union{Matrix,RkMatrix,HMatrix}) --> C

Two argument multiplication. When operating on `Union{Matrix,RkMatrix,HMatrix}`, the result `C` is returned in the *natural format*, as
described in the table below:

| `(*)(A,B)` | `B::Matrix`  | `B::RkMatrix`   | `B::HMatrix` |
|:-----|:---:|:---:|:---:|
|`A::Matrix`  | `C::Matrix`  | `C::RkMatrix` | `C::Matrix` |
|`A::RkMatrix`  | `C::RkMatrix` | `C::RkMatrix` | `C::RkMatrix` |
|`A::HMatrix`  | `C::Matrix` | `C::RkMatrix` | `C::HMatrix` |
"""
function (*)(M::Matrix,R::RkMatrix)
    tmp = M*R.A
    return RkMatrix(tmp,copy(R.B))
end

#1.3
# correctly falls back to mul!

#2.1
function (*)(R::RkMatrix,M::Matrix)
    tmp = adjoint(M)*R.B
    return RkMatrix(copy(R.A),tmp)
end

#2.2
function (*)(R::RkMatrix,S::RkMatrix)
    r,s = rank(R), rank(S)
    tmp = R.Bt*S.A
    return r<s ? RkMatrix(copy(R.A),S.B*adjoint(tmp)) : RkMatrix(R.A*tmp,copy(S.B))
end

#2.3
function (*)(R::RkMatrix,H::HMatrix)
    tmp = adjoint(H)*R.B
    return RkMatrix(copy(R.A),tmp)
end

#3.1
# correctly falls back to mul!

#3.2
function (*)(H::AbstractHMatrix,R::RkMatrix)
    tmp = H*R.A
    return RkMatrix(tmp,copy(R.B))
end

#3.3
function (*)(A::AbstractHMatrix,B::AbstractHMatrix)
    error("multiplying hierarchical matrix requires specifying a target hierarchical structure")
end

################################################################################
## mul!(C,A,B,a,b) :  C <-- b*C + a*A*B
## For the mul!(C,A,B,a,b) function, there are 3^3 = 27 cases to be considered
## depending on the types of C, A, and B. We will list these cases by by x.x.x
## where 1 means a full matrix, 2 a sparse matrix, and 3 a hierarhical matrix.
## E.g. case 1.2.1 means C is full, A is sparse, B is full
################################################################################

################################################################################
## 1.1.1
################################################################################
# default multiplication of dense matrices

################################################################################
## 1.1.2
################################################################################
function mul!(C::AbstractMatrix,M::AbstractMatrix,R::AbstractRkMatrix,a::Number,b::Number,buffer=similar(C,size(M,1),rank(R)))
    mul!(buffer,M,R.A)
    mul!(C,buffer,R.Bt,a,b)
    return C
end

################################################################################
## 1.1.3
################################################################################x
function mul!(C::AbstractMatrix,M::AbstractMatrix,H::AbstractHMatrix,a::Number,b::Number)
    rmul!(C,b)
    if hasdata(H)
        data = getdata(H)
        mul!(C,M,data,a,true)
    end
    for child in getchildren(H)
        shift  = pivot(H) .- 1
        irange = rowrange(child) .- shift[1]
        jrange = colrange(child) .- shift[2]
        Cview  = view(C,:,jrange)
        Mview  = view(M,:,irange)
        mul!(Cview,Mview,child,a,true)
    end
    return C
end

function _multiply_leaf!(C::Matrix,M::Matrix,H::AbstractHMatrix,a::Number,b::Number)
    irange = rowrange(H)
    jrange = colrange(H)
    Cview  = uview(C,:,jrange)
    Mview  = uview(M,:,irange)
    data   = getdata(H)
    mul!(Cview,Mview,data,a,b)
    return C
end

################################################################################
## 1.2.1
################################################################################
function mul!(C::AbstractMatrix,R::AbstractRkMatrix,M::AbstractMatrix,a::Number,b::Number)
    buffer = similar(C,rank(R),size(M,2))
    mul!(buffer,R.Bt,M)
    mul!(C,R.A,buffer,a,b)
end
function mul!(C::AbstractMatrix,adjR::Adjoint{<:Any,<:AbstractRkMatrix},M::AbstractMatrix,a::Number,b::Number)
    R = adjR.parent
    tmp = adjoint(R.A)*M
    mul!(C,R.B,tmp,a,b)
    return C
end

################################################################################
## 1.2.2
################################################################################
function mul!(C::AbstractMatrix,R::AbstractRkMatrix,S::AbstractRkMatrix,a::Number,b::Number)
    tmp = R*S
    mul!(C,tmp.A,tmp.Bt,a,b)
    return C
end

################################################################################
## 1.2.3
################################################################################
function mul!(C::AbstractMatrix,R::AbstractRkMatrix,H::AbstractHMatrix,a::Number,b::Number)
    tmp = R*H
    mul!(C,tmp.A,tmp.Bt,a,b)
    return C
end

################################################################################
## 1.3.1
################################################################################
function mul!(C::AbstractMatrix,H::AbstractHMatrix,M::AbstractMatrix,a::Number,b::Number)
    rmul!(C,b)
    if hasdata(H)
        data = getdata(H)
        mul!(C,data,M,a,true)
        # _multiply_leaf!(C,H,M,a,true)
    end
    for child in getchildren(H)
        shift  = pivot(H) .- 1
        irange = rowrange(child) .- shift[1]
        jrange = colrange(child) .- shift[2]
        Cview  = view(C,irange,:)
        Mview  = view(M,jrange,:)
        mul!(Cview,child,Mview,a,true)
    end
    return C
end

@inline function _multiply_leaf!(C::Matrix,H::AbstractHMatrix,M::Matrix,a,b)
    irange = rowrange(H)
    jrange = colrange(H)
    Cview  = view(C,irange,:)
    Mview  = view(M,jrange,:)
    data   = getdata(H)
    mul!(Cview,data,Mview,a,b)
    return C
end

function mul!(C::AbstractMatrix,adjH::Adjoint{<:Any,<:AbstractHMatrix},M::AbstractMatrix,a::Number,b::Number)
    rmul!(C,b)
    if hasdata(adjH)
        data = getdata(adjH)
        mul!(C,data,M,a,true)
    end
    for child in getchildren(adjH)
        shift  = pivot(adjH) .- 1
        irange = rowrange(child) .- shift[1]
        jrange = colrange(child) .- shift[2]
        Cview  = view(C,irange,:)
        Mview  = view(M,jrange,:)
        mul!(Cview,child,Mview,a,true)
    end
    return C
end

@inline function _multiply_leaf!(C::Matrix,adjH::Adjoint{<:Any,<:AbstractHMatrix},M::Matrix,a,b)
    irange = rowrange(adjH)
    jrange = colrange(adjH)
    Cview  = view(C,irange,:)
    Mview  = view(M,jrange,:)
    data   = getdata(adjH)
    mul!(Cview,data,Mview,a,b)
    return C
end

################################################################################
## 1.3.2
################################################################################
function mul!(C::AbstractMatrix,H::AbstractHMatrix,R::AbstractRkMatrix,a::Number,b::Number)
    buffer=similar(C,size(H,1),rank(R))
    mul!(buffer,H,R.A)
    mul!(C,buffer,R.Bt,a,b)
    return C
end

################################################################################
## 1.3.3 (should never arise in practice, thus sloppy implementation)
################################################################################
function mul!(C::AbstractMatrix,H::AbstractHMatrix,S::AbstractHMatrix,a::Number,b::Number)
    @debug "1.3.3: this case should not arise"
    mul!(C,H,Matrix(S),a,b)
    return C
end

################################################################################
## 2.1.1
################################################################################
# FIXME: should this case be considered?
function mul!(C::AbstractRkMatrix,M::AbstractMatrix,F::AbstractMatrix,a::Number,b::Number)
    @debug "2.1.1: this case should not arise"
    buffer=similar(C,size(M,1),size(F,2))
    mul!(buffer,M,F)
    axpby!(a,buffer,b,C)
    return C
end

################################################################################
## 2.1.2
################################################################################
function mul!(C::AbstractRkMatrix,M::AbstractMatrix,R::AbstractRkMatrix,a::Number,b::Number)
    tmp = M*R
    axpby!(a,tmp,b,C)
    return C
end

################################################################################
## 2.1.3
################################################################################
function mul!(C::AbstractRkMatrix,M::AbstractMatrix,H::AbstractHMatrix,a::Number,b::Number)
    @debug "2.1.3: this case should not arise"
    buffer = similar(C,size(M,1),size(H,2))
    mul!(buffer,M,H)
    axpby!(a,buffer,b,C)
    return C
end

################################################################################
## 2.2.1
################################################################################
function mul!(C::AbstractRkMatrix,R::AbstractRkMatrix,M::AbstractMatrix,a::Number,b::Number)
    tmp = R*M
    axpby!(a,tmp,b,C)
    return C
end

################################################################################
## 2.2.2
################################################################################
function mul!(C::AbstractRkMatrix,R::AbstractRkMatrix,S::AbstractRkMatrix,a::Number,b::Number)
    tmp = R*S
    axpby!(a,tmp,b,C)
    return C
end

################################################################################
## 2.2.3
################################################################################
function mul!(C::AbstractRkMatrix,R::AbstractRkMatrix,H::AbstractHMatrix,a::Number,b::Number)
    tmp = R*H
    axpby!(a,tmp,b,C)
    return C
end


################################################################################
## 2.3.1
################################################################################
function mul!(C::AbstractRkMatrix,H::AbstractHMatrix,M::AbstractMatrix,a::Number,b::Number)
    @debug "2.3.1: this case should not arise"
    T = promote_type(eltype(H),eltype(M))
    buffer = Matrix{T}(undef,size(H,1),size(M,2))
    mul!(buffer,H,M)
    axpby!(a,buffer,b,C)
    return C
end

################################################################################
## 2.3.2
################################################################################
function mul!(C::AbstractRkMatrix,H::AbstractHMatrix,R::AbstractRkMatrix,a::Number,b::Number)
    tmp = H*R
    axpby!(a,tmp,b,C)
    return C
end

################################################################################
## 2.3.3
################################################################################
function mul!(C::AbstractRkMatrix,A::AbstractHMatrix,B::AbstractHMatrix,a::Number,b::Number)
    # tp = TwoProd(A,B)
    # tmp = compress(tp,PartialACA(atol=1e-6))
    # R =
    # return axpby!(a,tmp,b,C)
    rmul!(C,b)
    if !isleaf(A) && !isleaf(B)
        m,n    = blocksize(A,1), blocksize(B,2)
        block  = Matrix{typeof(C)}(undef,m,n)
        for i = 1:m
            for j = 1:n
                p = size(getblock(A,i,1),1)
                q = size(getblock(B,1,j),2)
                block[i,j] = zero(typeof(C),p,q)
                for k = 1:blocksize(A,2)
                    mul!(block[i,j],getblock(A,i,k),getblock(B,k,j),true,true)
                end
            end
        end
        R = _gather(block)
        axpby!(a,R,true,C)
        # C.A = R.A
        # C.B = R.B
    else
        Adata = isleaf(A) ? A.data : A
        Bdata = isleaf(B) ? B.data : B
        mul!(C,Adata,Bdata,a,true)
    end
    return C
end

################################################################################
## 3.1.1
################################################################################
function mul!(C::AbstractHMatrix,M::AbstractMatrix,F::AbstractMatrix,a::Number,b::Number)
    tmp = M*F
    if hasdata(C)
        axpby!(a,tmp,b,C.data)
    else
        rmul!(C,b)
        C.data = rmul!(tmp,a)
    end
    return C
end

################################################################################
## 3.1.2
################################################################################
function mul!(C::AbstractHMatrix,M::AbstractMatrix,R::AbstractRkMatrix,a::Number,b::Number)
    tmp = M*R
    if hasdata(C)
        axpby!(a,tmp,b,C.data)
    else
        rmul!(C,b)
        C.data = rmul!(tmp,a)
    end
    return C
end

################################################################################
## 3.1.3
################################################################################
function mul!(C::AbstractHMatrix,M::AbstractMatrix,H::AbstractHMatrix,a::Number,b::Number)
    tmp = M*H
    if hasdata(C)
        axpby!(a,tmp,b,C.data)
    else
        rmul!(C,b)
        C.data = rmul!(tmp,a)
    end
    return C
end

################################################################################
## 3.2.1
################################################################################
function mul!(C::AbstractHMatrix,R::AbstractRkMatrix,M::AbstractMatrix,a::Number,b::Number)
    tmp = R*M
    if hasdata(C)
        axpby!(a,tmp,b,C.data)
    else
        rmul!(C,b)
        C.data = rmul!(tmp,a)
    end
    return C
end

################################################################################
## 3.2.2
################################################################################
function mul!(C::AbstractHMatrix,R::AbstractRkMatrix,M::AbstractRkMatrix,a::Number,b::Number)
    tmp = R*M
    if hasdata(C)
        axpby!(a,tmp,b,C.data)
    else
        rmul!(C,b)
        C.data = rmul!(tmp,a)
    end
    return C
end

################################################################################
## 3.2.3
################################################################################
function mul!(C::AbstractHMatrix,R::AbstractRkMatrix,H::AbstractHMatrix,a::Number,b::Number)
    tmp = R*H
    if hasdata(C)
        axpby!(a,tmp,b,C.data)
    else
        rmul!(C,b)
        C.data = rmul!(tmp,a)
    end
    return C
end

################################################################################
## 3.3.1
################################################################################
function mul!(C::AbstractHMatrix,H::AbstractHMatrix,M::Matrix,a::Number,b::Number)
    @debug "3.3.1: this case should not arise"
    tmp = H*M
    if hasdata(C)
        axpby!(a,tmp,b,C.data)
    else
        rmul!(C,b)
        C.data = rmul!(tmp,a)
    end
    return C
end

################################################################################
## 3.3.2
################################################################################
function mul!(C::AbstractHMatrix,H::AbstractHMatrix,R::AbstractRkMatrix,a::Number,b::Number)
    tmp = H*R
    if hasdata(C)
        axpby!(a,tmp,b,C.data)
    else
        rmul!(C,b)
        C.data = rmul!(tmp,a)
    end
    return C
end

################################################################################
## 3.3.3
################################################################################
function mul!(C::AbstractHMatrix,A::AbstractHMatrix,B::AbstractHMatrix,a::Number,b::Number)
    rmul!(C,b)
    if isleaf(A) || isleaf(B) || isleaf(B)
        _mul!(C,A,B,a,true)
    else
        # multiply the children
        ni,nj = blocksize(C)
        _ ,nk = blocksize(A)
        for i=1:ni
            for j=1:nj
                for k=1:nk
                    mul!(getblock(C,i,j),getblock(A,i,k),getblock(B,k,j),a,true)
                end
            end
        end
    end
    return C
end

# terminal case which dynamically dispatches to appropriate method
function _mul!(C::AbstractHMatrix,A::AbstractHMatrix,B::AbstractHMatrix,a::Number,b::Number)
    @assert isleaf(A) || isleaf(B) || isleaf(B)
    Cdata  = hasdata(C) ? getdata(C) : C
    Adata  = hasdata(A) ? getdata(A) : A
    Bdata  = hasdata(B) ? getdata(B) : B
    mul!(Cdata,Adata,Bdata,a,true)
    flush_tree!(C)
    return C
end

################################################################################
## FLUSH_TO_CHILDREN
################################################################################
function flush_to_children!(H::HMatrix)
    hasdata(H) && !isleaf(H) || (return H)
    R     = getdata(H)
    add_to_children!(H,R)
    setdata!(H,())
    return H
end

function add_to_children!(H,R::RkMatrix)
    shift = pivot(H) .- 1
    for block in getchildren(H)
        irange     = rowrange(block) .- shift[1]
        jrange     = colrange(block) .- shift[2]
        bdata      = getdata(block)
        tmp        = RkMatrix(R.A[irange,:],R.B[jrange,:])
        if bdata === ()
            setdata!(block,tmp)
        else
            axpby!(true,tmp,true,bdata)
        end
    end
end

function add_to_children!(H,M::Matrix)
    shift = pivot(H) .- 1
    for block in getchildren(H)
        irange     = rowrange(block) .- shift[1]
        jrange     = colrange(block) .- shift[2]
        bdata      = getdata(block)
        tmp        = M[irange,jrange]
        if bdata === ()
            setdata!(block,tmp)
        else
            axpby!(true,tmp,true,bdata)
        end
    end
end

################################################################################
## FLUSH_TREE
################################################################################
function flush_tree!(H::HMatrix)
    flush_to_children!(H)
    for block in getchildren(H)
        flush_tree!(block)
    end
    return H
end

################################################################################
## FLUSH_TO_LEAVES
################################################################################
function flush_to_leaves!(H::HMatrix)
    hasdata(H) && !isleaf(H) || (return H)
    R = getdata(H)
    add_to_leaves!(H,R)
    H.data = ()
    return H
end

function add_to_leaves!(H::HMatrix,R::RkMatrix)
    shift = pivot(H) .- 1
    for block in Leaves(H)
        irange     = rowrange(block) .- shift[1]
        jrange     = colrange(block) .- shift[2]
        bdata      = getdata(block)
        tmp        = RkMatrix(R.A[irange,:],R.B[jrange,:])
        if bdata === ()
            setdata!(block,tmp)
        else
            axpby!(true,tmp,true,bdata)
        end
    end
    return H
end

function add_to_leaves!(H::HMatrix,M::Matrix)
    shift = pivot(H) .- 1
    for block in Leaves(H)
        irange     = rowrange(block) .- shift[1]
        jrange     = colrange(block) .- shift[2]
        bdata      = getdata(block)
        tmp        = M[irange,jrange]
        if bdata === ()
            setdata!(block,tmp)
        else
            axpby!(true,tmp,true,bdata)
        end
    end
    return H
end
