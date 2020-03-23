function norm(R::RkFlexMatrix,p::Int=2)
    #TODO: improve this very rough estimation of the norm
    isempty(R) ? 0.0 : LinearAlgebra.norm(R.A[:,1],p)*LinearAlgebra.norm(R.B[:,1],p)
end
