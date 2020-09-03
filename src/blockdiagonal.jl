# https://bitbucket.org/cgeoga/blockdiagonal.jl/src/v0.0.1/src/BlockDiagonal.jl


struct BlockDiagonal{T, S} <: AbstractMatrix{T}
    blocks::Vector{S}
    partition::Vector{Int}
end

function BlockDiagonal(blocks)
    partition = Int[]
    push!(partition,0)
    for B in blocks
        push!(partition,last(partition)+size(B,1))
    end
    BlockDiagonal{eltype(first(blocks)), eltype(blocks)}(blocks, partition)
end

Base.size(B::BlockDiagonal) = (last(B.partition),last(B.partition))

function Base.getindex(B::BlockDiagonal,i::Int,j::Int)
    T  = eltype(B)
    p  = B.partition
    for k in 1:length(p)-1
        a,b = p[k]+1,p[k+1]
        if (a<=i<=b) && (a<=j<=b)
            M = B.blocks[k]
            iloc,jloc = i-a+1,j-a+1
            return M[iloc,jloc]
        end
    end
    return zero(T)
end

inv(B::BlockDiagonal)  = inv.(B.blocks) |> BlockDiagonal

function mul!(C::T,A::BlockDiagonal,B::T,a::Number,b::Number) where {T<:StridedVector}
    p = A.partition
    for k in 1:length(p)-1
        idxs = p[k]+1:p[k+1]
        M    = A.blocks[k]
        mul!(view(C,idxs),M,view(B,idxs),a,b)
    end
    return C
end

function mul!(C::T,A::BlockDiagonal,B::T,a::Number,b::Number) where {T<:StridedMatrix}
    p = A.partition
    for k in 1:length(p)-1
        idxs = p[k]+1:p[k+1]
        M    = A.blocks[k]
        mul!(view(C,idxs,:),M,view(B,idxs,:),a,b)
    end
    return C
end
