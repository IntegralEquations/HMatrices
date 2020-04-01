abstract type AbstractRkMatrix{T} <: AbstractMatrix{T} end

Base.size(rmat::AbstractRkMatrix)                                        = (size(rmat.A,1), size(rmat.B,1))
Base.isapprox(rmat::AbstractRkMatrix,B::AbstractArray,args...;kwargs...) = isapprox(Matrix(rmat),B,args...;kwargs...)
Base.getindex(rmat::AbstractRkMatrix,i::Int,j::Int)                      = sum(rmat.A[i,:].*rmat.Bt[:,j])

Base.getindex(rmat::AbstractRkMatrix,i::Int,::Colon) = getrow(rmat,i)
Base.getindex(rmat::AbstractRkMatrix,::Colon,j::Int) = getcol(rmat,j)
getcol(R,j)                                          = getcol!(Vector{eltype(R)}(undef,size(R,1)),R,j)
getrow(R,i)                                          = getrow!(Vector{eltype(R)}(undef,size(R,2)),R,i)
getcol!(col,R::AbstractRkMatrix,j::Int)              = mul!(col,R.A,uview(adjoint(R.B),:,j))
getrow!(row,R::AbstractMatrix,i)                     = mul!(row,R.B,uview(R.A,i,:))

function Base.getproperty(R::AbstractRkMatrix,s::Symbol)
    if  s == :Bt
        return adjoint(R.B)
    elseif  s == :At
        return adjoint(R.A)
    else
        return getfield(R,s)
    end
end

rank(M::AbstractRkMatrix) = size(M.A,2)

num_elements(R::AbstractRkMatrix)     = rank(R)*(sum(size(R)))
compression_rate(R::AbstractRkMatrix) = num_elements(R) / length(R)

"""
    RkFlexMatrix{T}

Similar to [`RkMatrix`](@ref), but internally stores the matrices `A` and `B` as `FlexMatrix{T}`.

See also: [`RkMatrix`](@ref), [`FlexMatrix`](@ref)
"""
struct RkFlexMatrix{T} <: AbstractRkMatrix{T}
    A::FlexMatrix{T}
    B::FlexMatrix{T}
    function RkFlexMatrix{T}(A::FlexMatrix{T},B::FlexMatrix{T}) where {T<:Number}
        @assert size(A,2) == size(B,2) "second dimension of `A` and `B` must match"
        m,r = size(A)
        n  = size(B,1)
        if  r*(m+n) >= m*n && m*n != 0
            @debug "Inefficient RkFlexMatrix: size(A)=$(size(A)), size(B)=$(size(B))"
        end
        new{T}(A,B)
    end
end
RkFlexMatrix(A::FlexMatrix{T},B::FlexMatrix{T}) where {T<:Number} = RkFlexMatrix{T}(A,B)
RkFlexMatrix(A::Vector{T},B::Vector{T}) where {T<:Vector}         = RkFlexMatrix(FlexMatrix(A),FlexMatrix(B))
RkFlexMatrix{T}(undef,m,n,r) where {T} = RkFlexMatrix(FlexMatrix{T}(undef,m,r),FlexMatrix{T}(undef,n,r))

function pushcross!(R::RkFlexMatrix,col,row)
    pushcol!(R.A,col)
    pushcol!(R.B,row)
    return R
end

"""
    RkMatrix{T}

Representation of a rank-`r` matrix ``M`` in the an outer product format:
```math M = AB^T ``` where ``A`` has size `m`×`r` and ``B`` has size `n`×`r`,
and ``B^T`` denotes the conjugate transpose (adjoint) of ``B``.
"""
mutable struct RkMatrix{T} <: AbstractRkMatrix{T}
    A::Matrix{T}
    B::Matrix{T}
    function RkMatrix{T}(A::Matrix,B::Matrix) where {T<:Number}
        @assert size(A,2) == size(B,2) "second dimension of `A` and `B` must match"
        m,r = size(A)
        n  = size(B,1)
        if  r*(m+n) >= m*n
            @debug "Inefficient RkMatrix: size(A)=$(size(A)), size(B)=$(size(B))"
        end
        new{T}(A,B)
    end
end
RkMatrix(A::Matrix{T},B::Matrix{T}) where {T} = RkMatrix{T}(A,B)
RkMatrix(A,B) = RkMatrix(promote(A,B)...)
RkMatrix(R::RkFlexMatrix) = RkMatrix(Matrix(R.A),Matrix(R.B))

function rkmatrix(F::LinearAlgebra.SVD)
    A  = F.U*LinearAlgebra.Diagonal(F.S)
    B  = copy(F.V)
    return RkMatrix(A,B)
end

function rkmatrix!(F::LinearAlgebra.SVD)
    A  = rmul!(F.U,LinearAlgebra.Diagonal(F.S))
    B  = F.V
    return RkMatrix(A,B)
end

function Base.hcat(M1::RkMatrix{T},M2::RkMatrix{T}) where {T}
    m,n  = size(M1)
    s,t  = size(M2)
    (m == s) || throw(ArgumentError("number of rows of each array must match: got  ($m,$s)"))
    r1   = size(M1.A,2)
    r2   = size(M2.A,2)
    A    = hcat(M1.A,M2.A)
    B1   = vcat(M1.B,zeros(T,t,r1))
    B2   = vcat(zeros(T,n,r2),M2.B)
    B    = hcat(B1,B2)
    return RkMatrix(A,B)
end

function Base.vcat(M1::RkMatrix{T},M2::RkMatrix{T}) where {T}
    m,n  = size(M1)
    s,t  = size(M2)
    n == t || throw(ArgumentError("number of columns of each array must match (got  ($n,$t))"))
    r1   = size(M1.A,2)
    r2   = size(M2.A,2)
    A1   = vcat(M1.A,zeros(T,s,r1))
    A2   = vcat(zeros(T,m,r2),M2.A)
    A    = hcat(A1,A2)
    B    = hcat(M1.B,M2.B)
    return RkMatrix(A,B)
end

function RkMatrix(M::Matrix)
    T = eltype(M)
    m,n = size(M)
    if m > n
        B = Matrix{T}(I,n,n)
        R = RkMatrix(M,B)
    else
        A  = Matrix{T}(I,m,m)
        B = adjoint(M) |> Matrix
        R = RkMatrix(A,B)
    end
    return R
end


Base.rand(::Type{RkMatrix{T}},m::Int,n::Int,r::Int) where {T} = RkMatrix(rand(T,m,r),rand(T,n,r)) #useful for testing
Base.rand(::Type{RkMatrix},m::Int,n::Int,r::Int) = rand(RkMatrix{Float64},m,n,r)

Base.copy(R::RkMatrix) = RkMatrix(copy(R.A),copy(R.B))

Matrix(R::RkMatrix) = R.A*R.Bt
