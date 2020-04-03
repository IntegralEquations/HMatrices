const HLU = LU{<:Any,<:AbstractHMatrix}

function Base.getproperty(LU::HLU,s::Symbol)
    T = eltype(LU)
    if s == :L
        H = LU.factors
        hasdata(H) ? (data = Matrix(UnitLowerTriangular(H.data))) : (data = ())
        if has_children(H)
            #diagonal blocks
            block11 = HLU(H[1,1]).L # recurses
            block22 = HLU(H[2,2]).L # recurses
            #zeroblock
            zeroData = zeros(RkMatrix{T},size(H[1,2],1),size(H[1,2],2),0.0)
            block12  = HMatrix{T}(H[1,2].irange,H[1,2].jrange,zeroData,(),())
            #
            block21  =  H[2,1]
            subblocks = (block11,block21,block12,block22)
            newH = HMatrix{T}(H.irange,H.jrange,data,(block11,block21,block12,block22),())
            for block in subblocks
                block.parent = newH
            end
            return newH
        else
            return HMatrix{T}(H.irange,H.jrange,data,(),())
        end
    elseif s == :U
        H = LU.hmat
        hasdata(H) ? (data = Matrix(UpperTriangular(H.data))) : (data = ())
        if has_children(H)
            #diagonal blocks
            block11 = HLU(H[1,1]).U # recurses
            block22 = HLU(H[2,2]).U # recurses
            #zeroblock
            zeroData = zeros(RkMatrix{T},size(H[2,1],1),size(H[2,1],2),0.0)
            block21  = HMatrix{T}(H[2,1].irange,H[2,1].jrange,zeroData,(),())
            #
            block12  =  H[1,2]
            subblocks = (block11,block21,block12,block22)
            newH = HMatrix{T}(H.irange,H.jrange,data,(block11,block21,block12,block22),())
            for block in subblocks
                block.parent = newH
            end
            return newH
        else
            return HMatrix{T}(H.irange,H.jrange,data,(),())
        end
    else
        return getfield(LU,s)
    end
end


function lu!(M::HMatrix,pivot=Val(false))
    if isleaf(M)
        data = getdata(M)
        @assert data isa Matrix
        lu!(data,pivot)
    else
        @assert !hasdata(M)
        M11 = getblock(M,1,1); M12 = getblock(M,1,2); M21 = getblock(M,2,1); M22 = getblock(M,2,2)
        lu!(M11,pivot)# M11 will now store L11 and U11
        ldiv!(M11.L,M12)     # L11*U12 = M12, M12 will store U12
        rdiv!(M21,M11.U)     # L21*U11 = M21, M21 will store L21
        mul!(M22,M21,M12,-1,1) # M22 = M22 - M21*M12
        lu!(M22,pivot)
    end
    isroot(M) ? (return HLU(M)) : return M
end
lu(M::HMatrix) = lu!(deepcopy(M))
