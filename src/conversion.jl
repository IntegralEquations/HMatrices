function RkMatrix(H::HMatrix)
    children = getchildren(H)
    if isleaf(H)
        hasdata(H) && (data = getdata(H))
        if data isa RkMatrix
            R = data
        else
            R = RkMatrix(data)
        end
    else
        B = [RkMatrix(child) for child in getchildren(H)]
        R = _gather(B)
    end
    return R
end

# FIXME: this is an inefficient way of gathering the block as it allocates too
# many temporary RkMatrices if size(block) != (2,2)
function _gather(B::Matrix{<:RkMatrix})
    cols = [_gather_row(B[i,:]) for i=1:size(B,1)]
    return _gather_col(cols)
end

function _gather_row(B::Vector{<:RkMatrix})
    R = B[1]
    for i=2:length(B)
        R = hcat(R,B[i])
    end
    return R
end

function _gather_col(B::Vector{<:RkMatrix})
    R = B[1]
    for i=2:length(B)
        R = vcat(R,B[i])
    end
    return R
end
