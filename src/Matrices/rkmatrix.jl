"""
    RkMatrix{T}

Representation of a rank-`r` matrix ``M`` in the an outer product format:
```math M = AB^T ``` where ``A`` has size `m`×`r` and ``B`` has size `n`×`r`,
and ``B^T`` denotes the conjugate transpose (adjoint) of ``B``.
"""
struct RkMatrix{T} <: AbstractMatrix{T}
    A::Matrix{T}
    B::Matrix{T}
    function RkMatrix{T}(A::Matrix,B::Matrix) where {T<:Number}
        @assert size(A,2) == size(B,2) "second dimension of `A` and `B` must match"
        m,r = size(A)
        n  = size(B,1)
        if  r*(m+n) >= m*n
            @debug "Inefficient RkMatrix: size(A)=$(size(A)), size(Bt)=$(size(Bt))"
        end
        new{T}(A,B)
    end
end
RkMatrix(A::Matrix{T},Bt::Matrix{T}) where {T} = RkMatrix{T}(A,Bt)
RkMatrix(A,B) = RkMatrix(promote(A,B)...)

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

function rkmatrix(A::Vector{V1},B::Vector{V2}) where {V1<: AbstractVector, V2 <: AbstractVector}
    T1   = eltype(V1)
    T2   = eltype(V2)
    T   = promote_type(T1,T2)
    r   = length(A)
    n   = length(first(A))
    m   = length(first(B))
    Ap  = Matrix{T}(undef,n,r)
    Bp = Matrix{T}(undef,m,r)
    for n=1:r
        copyto!(view(Ap,:,n),A[n])
        copyto!(view(Bp,:,n),B[n])
    end
    RkMatrix(Ap,Bp)
end

Base.size(rmat::RkMatrix)                                        = (size(rmat.A,1), size(rmat.B,1))
Base.size(rmat::RkMatrix,i)                                      = size(rmat)[i]
Base.length(rmat::RkMatrix)                                      = prod(size(rmat))
Base.isapprox(rmat::RkMatrix,B::AbstractArray,args...;kwargs...) = isapprox(Matrix(rmat),B,args...;kwargs...)
Base.getindex(rmat::RkMatrix,i::Int,j::Int)                      = sum(rmat.A[i,:].*rmat.Bt[:,j])
Base.getindex(rmat::RkMatrix,I,J)                                = RkMatrix(rmat.A[I,:],rmat.B[J,:])

function Base.getproperty(R::RkMatrix,s::Symbol)
    if  s == :Bt
        return adjoint(R.B)
    elseif  s == :At
        return adjoint(R.A)
    else
        return getfield(R,s)
    end
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

Base.rand(::Type{RkMatrix{T}},m::Int,n::Int,r::Int) where {T} = RkMatrix(rand(T,m,r),rand(T,n,r)) #useful for testing
Base.rand(::Type{RkMatrix},m::Int,n::Int,r::Int) = rand(RkMatrix{Float64},m,n,r)

num_elements(R::RkMatrix)        = size(R.A,2)*(sum(size(R)))
compression_rate(R::RkMatrix)    = num_elements(R) / prod(size(R))
