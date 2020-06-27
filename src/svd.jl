"""
    svd(R::RkMatrix)

Compute the singular value decomposition of an `RkMatrix` by first doing a `qr`
of `R.A` and `R.B` followed by an `svd` of ``R_A*(R_{B})^T``
"""
function svd(M::RkMatrix)
    r      = rank(M)
    # treat weird case where it would be most efficient to convert first to a full matrix
    r > min(size(M)...) && return svd(Matrix(M))
    # qr part
    QA, RA = qr(M.A)
    QB, RB = qr(M.B)
    # svd part
    F      = svd(RA*adjoint(RB))
    # build U  and Vt
    U      = QA*F.U
    Vt     = F.Vt*adjoint(QB)
    return SVD(U,F.S,Vt) #create the SVD structure
end

function svd!(M::RkMatrix)
    r      = rank(M)
    QA, RA = qr!(M.A)
    QB, RB = qr!(M.B)
    F      = svd!(RA*adjoint(RB))
    U      = QA*F.U
    Vt     = F.Vt*adjoint(QB)
    return SVD(U,F.S,Vt) #create the SVD structure
end
# TODO: the qrsvd function above should be optmized. A few notes: (a) by default
# Julia calls a blocked version of LAPACK's qr algorithm, but it appears the
# nonblocked version is faster (b) the blocksize keyword has an important effect
# on performance, and therefore could be tuned

# a posteriori truncation
function trunc_svd(A::AbstractMatrix{T},tol) where {T}
    F = svd(A)
    r = findlast(x -> x>tol, F.S)
    return SVD(F.U[:,1:r],F.S[1:r],F.V[:,1:r])
end
