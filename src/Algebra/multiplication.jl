hmul!(args...;kwargs...)                      = hmul!(CPU1(),args...;kwargs...)
hmul!(resource::CPU1,args...;kwargs...)       = _cpu1_hmul!(args...,kwargs...)
hmul!(resource::CPUThreads,args...;kwargs...) = _cputhreads_hmul!(args...,kwargs...)

LinearAlgebra.mul!(C::AbstractVector,H::AbstractHierarchicalMatrix,F::AbstractVector,a::Number,b::Number) = hmul!(C,H,F,a,b)
LinearAlgebra.mul!(r::AbstractResource,args...;kwargs...) = hmul!(r,args...;kwargs...)

function _cpu1_hmul!(C::AbstractVector,H::HMatrix,F::AbstractVector,a::Number,b::Number)
    b==1 || rmul!(C,b)
    if isleaf(H)
        data = getdata(H)
        mul!(C,data,F,a,1)
    else
        children = getchildren(H)
        m,n      = size(children)
        shift    = _idx_pivot(H) .- 1
        for i=1:m
            for j=1:n
                block  = children[i,j]
                irange = rowrange(block) .- shift[1]
                Cview  = view(C,irange)
                jrange = colrange(block) .- shift[2]
                Fview  = view(F,jrange)
                _cpu1_hmul!(Cview,block,Fview,a,1)
            end
        end
    end
    return C
end

function _cputhreads_hmul!(C::AbstractVector,H::HMatrix,F::AbstractVector,a::Number,b::Number)
    b==1 || rmul!(C,b)
    if isleaf(H)
        data = getdata(H)
        mul!(C,data,F,a,1)
    else
        children = getchildren(H)
        m,n      = size(children)
        shift    = _idx_pivot(H) .- 1
        t        = Matrix{Task}(undef,m,n)
        for i=1:m
            # TODO: choose between a ReentrantLock or a Semaphore here
            mutex = ReentrantLock()
            # mutex = Base.Semaphore(1)
            for j=1:n
                block  = children[i,j]
                irange = rowrange(block) .- shift[1]
                Cview  = view(C,irange)
                jrange = colrange(block) .- shift[2]
                Fview  = view(F,jrange)
                t[i,j] = @spawn begin
                    lock(mutex)
                    # Base.acquire(mutex)
                    _cputhreads_hmul!(Cview,block,Fview,a,1)
                    # Base.release(mutex)
                    unlock(mutex)
                end
            end
        end
        waitall(t)
    end
    return C
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
