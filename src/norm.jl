# function norm(R::RkFlexMatrix,p::Int=2)
#     #FIXME: improve this very rough estimation of the norm
#     isempty(R) ? 0.0 : norm(R.A[:,1],p)*norm(R.B[:,1],p)
# end
