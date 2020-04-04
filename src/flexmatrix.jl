"""
    FlexMatrix{T}

Represents a matrix as a vector of vectors.
"""
struct FlexMatrix{T} <: AbstractMatrix{T}
    data::Vector{Vector{T}}
    function FlexMatrix{T}(data::Vector{Vector{T}}) where {T<:Number}
        @assert _ismatrix(data) "data must correspond to a matrix"
        new{T}(data)
    end
end
FlexMatrix(M::Vector{Vector{T}}) where {T} = FlexMatrix{T}(M)

FlexMatrix{T}(undef,m,n) where {T} = FlexMatrix([Vector{T}(undef,m) for _=1:n])

function FlexMatrix(M::Matrix)
    m,n  = size(M)
    T    = eltype(M)
    F    = FlexMatrix{T}(undef,m,n)
    for n=1:size(F,2)
        copyto!(F[:,n],view(M,:,n))
    end
    return F
end

Base.isempty(M::FlexMatrix)    = isempty(M.data)
Base.size(M::FlexMatrix)       = isempty(M) ? (0,0) : (length(first(M.data)), length(M.data))
Base.getindex(A::FlexMatrix,i::Int,j::Int)  = A.data[j][i]
Base.getindex(A::FlexMatrix,::Colon,j::Int) = A.data[j]

function Base.Matrix(F::FlexMatrix)
    T = eltype(F)
    M = Matrix{T}(undef,size(F))
    for n=1:size(M,2)
        copyto!(view(M,:,n),F[:,n])
    end
    return M
end

function pushcol!(F::FlexMatrix,vs...)
    # TODO: assert hat the vectors vs being added are of the correct length
    push!(F.data,vs...)
    return F
end

function _ismatrix(A::Vector{<:Vector})
    isempty(A) && return true
    n = length(first(A))
    for v in A
        length(v) == n || (return false)
    end
    return true
end
