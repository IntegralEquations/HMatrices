################################################################################
## INV
################################################################################
function inv!(M::HMatrix,X::HMatrix=zero(M))
    if isleaf(M)
        data = getdata(M)
        setdata!(M,inv(data))
    else
        @assert !hasdata(M) #only leaves are allowed to have data for the inversion
        #recursion
        inv!(getblock(M,1,1),getblock(X,1,1))
        #update
        mul!(getblock(X,1,2),getblock(M,1,1),getblock(M,1,2),-1,false)
        flush_tree!(getblock(X,1,2))

        mul!(getblock(X,2,1),getblock(M,2,1),getblock(M,1,1),true,false)
        flush_tree!(getblock(X,2,1))

        mul!(getblock(M,2,2),getblock(M,2,1),getblock(X,1,2),true,true)
        flush_tree!(getblock(M,2,2))

        #recursion
        inv!(getblock(M,2,2),getblock(X,2,2))

        #update
        mul!(getblock(M,1,2),getblock(X,1,2),getblock(M,2,2),1,0)
        flush_tree!(getblock(M,1,2))

        mul!(getblock(M,1,1),getblock(M,1,2),getblock(X,2,1),-1,true)
        flush_tree!(getblock(M,1,1))

        mul!(getblock(M,2,1),getblock(M,2,2),getblock(X,2,1),-1,false)
        flush_tree!(getblock(M,2,1))
    end
    return M
end

inv(M::HMatrix) = inv!(deepcopy(M))
