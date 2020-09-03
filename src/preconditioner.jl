# abstract type AbstractPreconditioner{T} <: AbstractMatrix{T} end
abstract type AbstractPreconditioner{T} end

Base.size(P::AbstractPreconditioner,args...) = size(P.matrix,args...)
Base.getindex(P::AbstractPreconditioner,args...) = getindex(P.matrix,args...)

struct DiagonalPreconditioner{T} <: AbstractPreconditioner{T}
    matrix::Diagonal{T,Vector{T}}
    inverse::Diagonal{T,Vector{T}}
end

DiagonalPreconditioner(D::Diagonal) = DiagonalPreconditioner(D,inv(D))
DiagonalPreconditioner(v::Vector)   = DiagonalPreconditioner(Diagona(v))
DiagonalPreconditioner(H::AbstractHMatrix)  = DiagonalPreconditioner(Diagonal(H))

struct BlockDiagonalPreconditioner{T} <: AbstractPreconditioner{T}
    matrix::BlockDiagonal{T,Matrix{T}}
    inverse::BlockDiagonal{T,Matrix{T}}
    lu::Vector{LU{T,Matrix{T}}}
end

BlockDiagonalPreconditioner(B::BlockDiagonal)         = BlockDiagonalPreconditioner(B,inv(B),lu.(B.blocks))
BlockDiagonalPreconditioner(blocks::Vector{<:Matrix}) = BlockDiagonalPreconditioner(BlockDiagonal(blocks))
BlockDiagonalPreconditioner(H::AbstractHMatrix)       = BlockDiagonalPreconditioner(BlockDiagonal(H))

ldiv!(y,P::AbstractPreconditioner,x) = mul!(y,P.inverse,x)

ldiv!(P::AbstractPreconditioner,x) = mul!(x,P.inverse,deepcopy(x))

function ldiv!(P::DiagonalPreconditioner,x)
    x .= P.inverse.diag .* x
    return x
end

function ldiv!(P::BlockDiagonalPreconditioner,x)
    A = P.matrix
    p = A.partition
    for k in 1:length(p)-1
        idxs = p[k]+1:p[k+1]
        F    = P.lu[k]
        ldiv!(F,view(x,idxs))
    end
    return x
end
