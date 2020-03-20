hmul!(args...;kwargs...)                      = hmul!(CPU1(),args...;kwargs...)
hmul!(resource::CPU1,args...;kwargs...)       = _cpu1_hmul!(args...,kwargs...)
hmul!(resource::CPUThreads,args...;kwargs...) = _cputhreads_hmul!(args...,kwargs...)

LinearAlgebra.mul!(C::AbstractVector,H::AbstractHierarchicalMatrix,F::AbstractVector,a::Number,b::Number) = hmul!(C,H,F,a,b)
LinearAlgebra.mul!(r::AbstractResource,args...;kwargs...) = hmul!(r,args...;kwargs...)

function _cpu1_hmul!(C::AbstractVector,H::HMatrix,F::AbstractVector,a::Number=true,b::Number=false)
    b==1 || rmul!(C,b)
    if isleaf(H)
        _multiply_leaf!(C,H,F,a,1)
    else
        for child in getchildren(H)
            _cpu1_hmul!(C,child,F,a,1)
        end
    end
    return C
end

# function _cputhreads_hmul!(C::AbstractVector,H::HMatrix,F::AbstractVector,a::Number=true,b::Number=false)
#     b==1 || rmul!(C,b)
#     if isleaf(H)
#         _multiply_leaf!(C,H,F,a,1)
#     else
#         children = getchildren(H)
#         m,n      = size(children)
#         t        = Matrix{Task}(undef,m,n)
#         for i=1:m
#             # NOTE: either a ReentrantLock or a Semaphore can be used here. They
#             # seem to perform about the same in this context
#             mutex = ReentrantLock()
#             for j=1:n
#                 block  = children[i,j]
#                 t[i,j] = @spawn begin
#                     lock(mutex)
#                     _cputhreads_hmul!(C,block,F,a,1)
#                     unlock(mutex)
#                 end
#             end
#         end
#         waitall(t)
#     end
#     return C
# end

function _cputhreads_hmul!(C::AbstractVector,H::T,F::AbstractVector,a::Number=true,b::Number=false) where {T<:HMatrix}
    rmul!(C,b)
    nthreads = Threads.nthreads()
    leaves   = LeavesChannel(H,Inf)
    sort!(leaves.data,lt = (a,b)->length(a)<length(b))
    y        = [zero(C) for _ = 1:nthreads]
    @sync for i in 1:nthreads
        @spawn begin
            _multiply_leaves!(y[i],leaves,F,a,1)
        end
    end
    #reduction stage
    for i=1:nthreads
        axpy!(1,y[i],C)
    end
    return C
end

function _multiply_leaves!(y,leaves,x,a,b)
    for leaf in leaves
        _multiply_leaf!(y,leaf,x,a,b)
    end
    return y
end

function _multiply_leaf!(C,H,F,a,b)
    irange = rowrange(H)
    jrange = colrange(H)
    Cview  = view(C,irange)
    Fview  = view(F,jrange)
    data   = getdata(H)
    mul!(Cview,data,Fview,a,b)
    return Cview
end

waitall(t) = [wait(ti) for ti in t]

function LinearAlgebra.mul!(C::AbstractVector,Rk::RkFlexMatrix,F::AbstractVector,a::Number,b::Number,buffer=similar(C,rank(Rk)))
    buffer = mul!(buffer,Rk.Bt,F)
    mul!(C,Rk.A,buffer,a,b)
end

# memory_overlap(A,R::RkFlexMatrix) = false
# memory_overlap(R::RkFlexMatrix,A) = false
# memory_overlap(R::RkFlexMatrix,A::RkFlexMatrix) = false
# memory_overlap(R::FlexMatrix,A)   = false
# memory_overlap(A,R::FlexMatrix)   = false
# function plan_mul!(C::AbstractVector,H::HMatrix,F::AbstractVector,a::Number,b::Number,tg::TaskGraph=TaskGraph())
#     # @addop tg LinearAlgebra.rmul!(C,b)
#     cl = @codelet LinearAlgebra.rmul!(C,b)
#     HScheduler.addvertex!(tg,cl)
#     shift = _idx_pivot(H) .- 1
#     for block in AbstractTrees.Leaves(H)

#         irange = rowrange(block) .- shift[1]
#         jrange = colrange(block) .- shift[2]
#         data = getdata(block)
#         Cview = view(C,irange)
#         Fview = view(F,jrange)
#         cl = @codelet LinearAlgebra.mul!(Cview,data,Fview,a,1)
#         HScheduler.addvertex!(tg,cl)
#     end
#     cl = Codelet(func=(C)->(C),data=(C,),access_mode=(HScheduler.READ,))
#     HScheduler.addvertex!(tg,cl)
#     return tg
# end
