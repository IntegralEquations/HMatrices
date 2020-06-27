abstract type AbstractCompressor end

"""
    ACA <: AbstractCompressor

Adaptive cross approximation compressor with full pivoting for finding cross. Requires evaluation of entire matrix, but
is guaranteed  to work.
"""
@Base.kwdef struct ACA <: AbstractCompressor
    atol::Float64 = 0
    rank::Int     = typemax(Int)
    rtol::Float64 = atol>0 || rank<typemax(Int) ? 0 : sqrt(eps(Float64))
    p::Int = 2
end

"""
    PartialACA <: AbstractCompressor

Adaptive cross approximation compressor with partial pivoting for finding cross. Does not require evaluation of entire matrix, but is not guaranteed to converge either.
"""
@Base.kwdef struct PartialACA <: AbstractCompressor
    atol::Float64 = 0
    rank::Int     = typemax(Int)
    rtol::Float64 = atol>0 || rank<typemax(Int) ? 0 : sqrt(eps(Float64))
end

"""
    TSVD

Truncated svd
"""
@Base.kwdef struct TSVD
    atol::Float64 = 0
    rank::Int     = typemax(Int)
    rtol::Float64 = atol>0 || rank<typemax(Int) ? 0 : sqrt(eps(Float64))
end


function (aca::ACA)(K,irange::UnitRange,jrange::UnitRange)
    M  = K[irange,jrange] #computes the entire matrix.
    _aca_full!(M,aca.atol,aca.rank,aca.rtol,x->norm(x,aca.p))
end

function _aca_full!(M, atol, rmax, rtol, norm)
    T   = eltype(M)
    m,n = size(M)
    R   = RkFlexMatrix{T}(undef,m,n,0)
    er = Inf
    exact_norm = norm(M) #exact norm
    while er > max(atol,rtol*exact_norm) && rank(R) < rmax
        (i,j) = argmax(abs.(M)).I
        δ       = M[i,j]
        if δ == 0
            return R
        else
            a = M[:,j]
            b = conj(M[i,:])
            rdiv!(a,δ)
            pushcross!(R,a,b)
            axpy!(-1,a*adjoint(b),M) # M <-- M - col*row'
            er = norm(M) # exact error
        end
    end
    return RkMatrix(R)
end

function (paca::PartialACA)(K,irange::UnitRange,jrange::UnitRange)
    _aca_partial(K,irange,jrange,paca.atol,paca.rank,paca.rtol)
end

function _aca_partial(K,irange,jrange,atol,rmax,rtol)
    ishift,jshift = irange.start-1, jrange.start-1 #maps global indices to local indices
    T   = Base.eltype(K)
    m,n = length(irange),length(jrange)
    rmax = min(rmax,min(m,n))
    R   = RkFlexMatrix{T}(undef,m,n,0)
    A   = R.A
    B   = R.B
    I   = [true for i = 1:m]
    J   = [true for i = 1:n]
    i   = 1
    er  = Inf
    est_norm = 0 #approximate norm of K
    while er > max(atol,rtol*est_norm) && rank(R) < rmax
        # remove index i from allowed row
        I[i] = false
        # compute next row by b <-- conj(K[i+ishift,jrange] - R[i,:])
        b    = conj!(K[i+ishift,jrange])
        for k = 1:rank(R)
            axpy!(-conj(A[i,k]),B[:,k],b)
        end
        j    = _nextcol(b,J)
        δ    = b[j]
        if δ == 0
            i = findfirst(x->x==true,J)
            isnothing(i) && error("aca did not converge")
        else # δ != 0
            rdiv!(b,δ) # b <-- b/δ
            J[j] = false
            # compute next col by a <-- K[irange,j+jshift] - R[:,j]
            a    = K[irange,j+jshift]
            for k = 1:rank(R)
                axpy!(-conj(B[j,k]),A[:,k],a)
            end
            pushcross!(R,a,b)
            er       = norm(a)*norm(b) # estimate the error by || R_{k} - R_{k-1} || = ||a|| ||b||
            est_norm = _update_frob_norm(est_norm,R) # estimate the norm by || K || ≈ || R_k ||
            i        = _nextrow(a,I)
        end
    end
    return RkMatrix(R)
end

@inline function _update_frob_norm(cur,R)
    k = rank(R)
    A,B = R.A, R.B
    a = A[:,end]
    b = B[:,end]
    out = norm(a)^2 * norm(b)^2
    for l=1:k-1
        out += 2*real(dot(A[:,l],a)*conj(dot(B[:,l],b)))
    end
    return sqrt(cur^2 + out)
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


################################################################################
## Truncated SVD
###############################################################################
compress(R::RkMatrix,tsvd::TSVD) = compress!(deepcopy(R),tsvd)
function compress!(R::RkMatrix,tsvd::TSVD)
    m,n = size(R)
    F   = svd!(R)
    enorm = F.S[1]
    r = findlast(x -> x>max(tsvd.atol,tsvd.rtol*enorm), F.S)
    r = min(r,tsvd.rank)
    if m<n
        R.A = F.U[:,1:r]*Diagonal(F.S[1:r])
        R.B = F.V[:,1:r]
    else
        R.A = F.U[:,1:r]
        R.B = F.V[:,1:r]*adjoint(Diagonal(F.S[1:r]))
    end
    return R
end

compress(M::Matrix,tsvd::TSVD) = compress!(copy(M),tsvd)
function compress!(M::Matrix,tsvd::TSVD)
    F = svd!(M)
    enorm = F.S[1]
    r = findlast(x -> x>max(tsvd.atol,tsvd.rtol*enorm), F.S)
    r = min(r,tsvd.rank)
    m,n = size(M)
    if m<n
        A = @views F.U[:,1:r]*Diagonal(F.S[1:r])
        B = F.V[:,1:r]
    else
        A = F.U[:,1:r]
        B = @views F.V[:,1:r]*adjoint(Diagonal(F.S[1:r]))
    end
    return RkMatrix(A,B)
end

function compress!(H::HMatrix,tsvd::TSVD)
    if isadmissible(H)
        compress!(getdata(H),tsvd)
    else
        for child in getchildren(H)
            compress(child,tsvd)
        end
    end
    return H
end
