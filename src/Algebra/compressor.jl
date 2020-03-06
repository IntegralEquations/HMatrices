abstract type AbstractCompressor end

compress(K,block,compressor::AbstractCompressor) = compress(K,rowrange(block),colrange(block),compressor)

"""
    ACA <: AbstractCompressor

Adaptive cross approximation compressor with full pivoting for finding cross. Requires evaluation of entire matrix, but
is guaranteed  to work.
"""
@Base.kwdef struct ACA
    atol::Float64 = 0
    rank::Int     = typemax(Int)
    rtol::Float64 = atol>0 || rank<typemax(Int) ? 0 : sqrt(eps(Float64))
    p::Int = 2
end

function compress(K,irange::UnitRange,jrange::UnitRange,aca::ACA)
    M  = K[irange,jrange] #computes the entire matrix.
    _aca_full!(M,aca.atol,aca.rank,aca.rtol,x->LinearAlgebra.norm(x,aca.p))
end

function _aca_full!(M, atol, rmax, rtol, norm)
    T   = eltype(M)
    m,n = size(M)
    R   = RkFlexMatrix{T}(undef,m,n,0)
    er = Inf
    exact_norm = norm(M) #exact norm
    while er > max(atol,rtol*exact_norm) && rank(R) < rmax
        @debug er
        (i,j) = argmax(abs.(M)).I
        δ       = M[i,j]
        if δ == 0
            return R
        else
            a = M[:,j]
            b = conj(M[i,:])
            LinearAlgebra.rdiv!(a,δ)
            pushcross!(R,a,b)
            LinearAlgebra.axpy!(-1,a*adjoint(b),M) # M <-- M - col*row'
            er = norm(M) # exact error
        end
    end
    return R
end

"""
    PartialACA <: AbstractCompressor

Adaptive cross approximation compressor with partial pivoting for finding cross. Does not require evaluation of entire matrix, but is not guaranteed to converge either.
"""
@Base.kwdef struct PartialACA
    atol::Float64 = 0
    rank::Int     = typemax(Int)
    rtol::Float64 = atol>0 || rank<typemax(Int) ? 0 : sqrt(eps(Float64))
    p::Int = 2
end

function compress(K,irange::UnitRange,jrange::UnitRange,aca::PartialACA)
    _aca_partial(K,irange,jrange,aca.atol,aca.rank,aca.rtol,x->LinearAlgebra.norm(x,aca.p))
end

function _aca_partial(K,irange,jrange,atol,rmax,rtol,norm)
    T   = Base.eltype(K)
    m,n = length(irange),length(jrange)
    R   = RkFlexMatrix{T}(undef,m,n,0)
    I   = [true for i = 1:m]
    J   = [true for i = 1:n]
    i   = 1
    er  = Inf
    est_norm = 0 #approximate norm of K
    while er > max(atol,rtol*est_norm) && rank(R) < rmax
        I[i] = false  # remove index i from allowed row
        b    = isempty(R) ? conj(K[i,jrange]) : conj(K[i,:] - R[i,:])
        j    = _nextcol(b,J)
        δ    = b[j]
        if δ == 0
            i = findfirst(x->x==true,J)
            isnothing(i) && error("aca did not converge")
        else # δ != 0
            LinearAlgebra.rdiv!(b,δ) # b <-- b/δ
            J[j] = false
            a    = isempty(R) ? K[:,j] : K[:,j] - R[:,j] # compute a column
            pushcross!(R,a,b)
            er       = norm(a)*norm(b) # approximate error
            est_norm = norm(R)
            i        = _nextrow(a,I)
        end
    end
    return R
end

function _nextcol(col,J)
    out = -1
    val = -Inf
    for n in 1:length(J)
        J[n] || continue
        tmp = abs(col[n])
        tmp < val && continue
        out = n
        val = tmp
    end
    return out
end
_nextrow(row,I) = _nextcol(row,I)

function LinearAlgebra.norm(R::RkFlexMatrix,p::Int=2)
    #TODO: improve this very rough estimation of the norm
    isempty(R) ? 0.0 : LinearAlgebra.norm(R.A[:,1],p)*LinearAlgebra.norm(R.B[:,1],p)
end
