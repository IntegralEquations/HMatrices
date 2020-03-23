################################################################################
## Utilities for working with trees
################################################################################
function depth(node,acc=0)
    if isroot(node)
        return acc
    else
        depth(getparent(node),acc+1)
    end
end

children(H::AbstractHierarchicalMatrix) = getchildren(H)
Base.eltype(::Type{<:TreeIterator{T}}) where {T<:AbstractHierarchicalMatrix}         = T
Base.IteratorEltype(::Type{<:TreeIterator{T}}) where {T<:AbstractHierarchicalMatrix} = Base.HasEltype()

