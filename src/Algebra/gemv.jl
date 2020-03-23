const RkType = Union{RkFlexMatrix,RkMatrix}

function mul!(C::AbstractVector,Rk::RkType,F::AbstractVector,a::Number,b::Number,buffer=similar(C,rank(Rk)))
    buffer = mul!(buffer,Rk.Bt,F)
    mul!(C,Rk.A,buffer,a,b)
end

# CPU1
function mul!(resources::CPU1,C::AbstractVector,H::AbstractHierarchicalMatrix,F::AbstractVector,a::Number,b::Number)
    b==1 || rmul!(C,b)
    for leaf in Leaves(H)
        _multiply_leaf!(C,leaf,F,a,1)
    end
    return C
end

# CPUThreads
function mul!(resources::CPUThreads,C::AbstractVector,H::AbstractHierarchicalMatrix,F::AbstractVector,a::Number=true,b::Number=false)
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
mul!(C::AbstractVector,H::AbstractHierarchicalMatrix,F::AbstractVector,a::Number,b::Number) = mul!(CPU1(),C,H,F,a,b)

# 
function _multiply_leaf!(C,H,F,a,b)
    irange = rowrange(H)
    jrange = colrange(H)
    Cview  = view(C,irange)
    Fview  = view(F,jrange)
    data   = getdata(H)
    mul!(Cview,data,Fview,a,b)
    return Cview
end

# NOTE: the commented code below avoids making copies at the cost of locking mutex to avoid race conditions. 
# TODO: decide on the best one.
# function mul!(resources::CPUThreads,C::AbstractVector,H::AbstractHierarchicalMatrix,F::AbstractVector,a::Number=true,b::Number=false)
#     b==1 || rmul!(C,b)
#     if isleaf(H)
#         _multiply_leaf!(C,H,F,a,1)
#     else
#         children = getchildren(H)
#         m,n      = size(children)
#         t        = Matrix{Task}(undef,m,n)
#         for i=1:m
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


