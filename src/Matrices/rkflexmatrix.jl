"""
    RkMatrixFlex{T}

Similar to [`RkMatrix`](@ref), but internally stores the matrices `A` and `B` as `FlexMatrix{T}`.

See also: [`RkMatrix`](@ref), [`FlexMatrix`](@ref)
"""
struct RkFlexMatrix{T} <: AbstractMatrix{T}
    A::FlexMatrix{T}
    B::FlexMatrix{T}
    function RkFlexMatrix{T}(A::FlexMatrix{T},B::FlexMatrix{T}) where {T<:Number}
        @assert size(A,2) == size(B,2) "second dimension of `A` and `B` must match"
        m,r = size(A)
        n  = size(B,1)
        if  r*(m+n) >= m*n
            @debug "Inefficient RkFlexMatrix: size(A)=$(size(A)), size(B)=$(size(B))"
        end
        new{T}(A,B)
    end
end
RkFlexMatrix(A::FlexMatrix{T},B::FlexMatrix{T}) where {T<:Number} = RkFlexMatrix{T}(A,B)

RkFlexMatrix{T}(undef,m,n,r) where {T} = RkFlexMatrix(FlexMatrix{T}(undef,m,r),FlexMatrix{T}(undef,n,r))

Base.size(rmat::RkFlexMatrix)                                        = (size(rmat.A,1), size(rmat.B,1))
Base.isapprox(rmat::RkFlexMatrix,B::AbstractArray,args...;kwargs...) = isapprox(Matrix(rmat),B,args...;kwargs...)
Base.getindex(rmat::RkFlexMatrix,i::Int,j::Int)                      = sum(rmat.A[i,:].*rmat.Bt[:,j])
Base.getindex(rmat::RkFlexMatrix,I,J)                                = RkFlexMatrix(rmat.A[I,:],rmat.B[J,:])

function Base.getproperty(R::RkFlexMatrix,s::Symbol)
    if  s == :Bt
        return adjoint(R.B)
    elseif  s == :At
        return adjoint(R.A)
    else
        return getfield(R,s)
    end
end

function pushcross!(R::RkFlexMatrix,col,row)
    @debug length(col), length(row), size(R)
    pushcol!(R.A,col)
    pushcol!(R.B,row)
    return R
end
