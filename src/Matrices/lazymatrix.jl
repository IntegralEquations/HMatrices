struct LazyMatrix{Tf,Tx,Ty,T} <: AbstractMatrix{T}
    f::Tf
    X::Tx
    Y::Ty
end
Base.size(lz::LazyMatrix)                   = length(lz.X), length(lz.Y)
Base.getindex(lz::LazyMatrix,i::Int,j::Int) =  lz.f(lz.X[i],lz.Y[j])

Base.conj(lz::LazyMatrix) = LazyMatrix(conj(lz.f),lz.X,lz.Y)

function LazyMatrix(f,X,Y)
    T = Base.promote_op(f,eltype(X),eltype(Y))
    LazyMatrix{typeof(f),typeof(X),typeof(Y),T}(f,X,Y)
end

#TODO: support algebra of LazyMatrix?
