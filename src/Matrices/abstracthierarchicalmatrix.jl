abstract type  AbstractHierarchicalMatrix{T} <: AbstractMatrix{T} end

function Base.getindex(H::AbstractHierarchicalMatrix,i::Int,j::Int)
    (i ∈ rowrange(H)) && (j ∈ colrange(H)) || throw(BoundsError(H,(i,j)))
    if isleaf(H)
        hasdata(H) || (return nothing)
        il,jl = _idx_global_to_local(i,j,H)
        return getdata(H)[il,jl]
    else
        for child in getchildren(H)
            if (i ∈ rowrange(child)) && (j ∈ colrange(child))
                return getindex(child,i,j)
            end
        end
    end
end

Base.size(H::AbstractHierarchicalMatrix) = length(rowrange(H)), length(colrange(H))

rowrange(H::AbstractHierarchicalMatrix) = H.rowrange
colrange(H::AbstractHierarchicalMatrix) = H.colrange
getchildren(H::AbstractHierarchicalMatrix)                 = H.children
setchildren!(H::AbstractHierarchicalMatrix,chd)            = (H.children = chd)
getparent(H::AbstractHierarchicalMatrix)        = H.parent
setparent!(H::AbstractHierarchicalMatrix,par)   = (H.parent = par)
getdata(H::AbstractHierarchicalMatrix)          = H.data
setdata!(H::AbstractHierarchicalMatrix,data)    = (H.data = data)

_idx_pivot(H::AbstractHierarchicalMatrix)   = rowrange(H).start, colrange(H).start
_idx_global_to_local(I,J,H::AbstractHierarchicalMatrix) = (I,J) .- _idx_pivot(H) .+ 1
isleaf(H::AbstractHierarchicalMatrix)       = getchildren(H) === ()
isroot(H::AbstractHierarchicalMatrix)       = getparent(H) === ()
hasdata(H::AbstractHierarchicalMatrix)      = getdata(H) !== ()

AbstractTrees.children(H::AbstractHierarchicalMatrix) = getchildren(H)
