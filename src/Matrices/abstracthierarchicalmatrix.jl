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

rowrange(H::AbstractHierarchicalMatrix)         = H.rowrange
colrange(H::AbstractHierarchicalMatrix)         = H.colrange
pivot(H::AbstractHierarchicalMatrix)            = (rowrange(H).start,colrange(H).start)
getchildren(H::AbstractHierarchicalMatrix)      = H.children
setchildren!(H::AbstractHierarchicalMatrix,chd) = (H.children = chd)
getparent(H::AbstractHierarchicalMatrix)        = H.parent
setparent!(H::AbstractHierarchicalMatrix,par)   = (H.parent   = par)
getdata(H::AbstractHierarchicalMatrix)          = H.data
setdata!(H::AbstractHierarchicalMatrix,data)    = (H.data     = data)
isadmissible(H::AbstractHierarchicalMatrix)     = H.admissible


_idx_pivot(H::AbstractHierarchicalMatrix)               = rowrange(H).start, colrange(H).start
_idx_global_to_local(I,J,H::AbstractHierarchicalMatrix) = (I,J) .- _idx_pivot(H) .+ 1
isleaf(H::AbstractHierarchicalMatrix)                   = getchildren(H) === ()
isroot(H::AbstractHierarchicalMatrix)                   = getparent(H) === ()
hasdata(H::AbstractHierarchicalMatrix)                  = getdata(H) !== ()

function compression_rate(H::AbstractHierarchicalMatrix)
    c = 0 # stored entries
    for leaf in Leaves(H)
        rel_size = length(leaf)/length(H)
        if isadmissible(leaf)
            c += compression_rate(leaf.data)*rel_size
        else
            c += rel_size
        end
    end
    return c
end

# Interface to AbstractTrees
children(H::AbstractHierarchicalMatrix) = getchildren(H)
Base.eltype(::Type{<:TreeIterator{T}}) where {T<:AbstractHierarchicalMatrix}         = T
Base.IteratorEltype(::Type{<:TreeIterator{T}}) where {T<:AbstractHierarchicalMatrix} = Base.HasEltype()

################################################################################
## Plot recipes
################################################################################
@recipe function f(hmat::HMatrix)
    legend --> false
    grid   --> false
    # aspect_ratio --> :equal
    yflip  := true
    seriestype := :shape
    linecolor  --> :black
    title  := "hmatrix"
    # all leaves
    for block in Leaves(hmat)
        @series begin
            if isadmissible(block)
                fillcolor  := :blue
                opacity    := compression_rate(getdata(block))
            else
                fillcolor  := :red
                opacity     = 0.3
            end
            pt1 = pivot(block)
            pt2 = pt1 .+ size(block) .- 1
            y1, y2 = pt1[1],pt2[1]
            x1, x2 = pt1[2],pt2[2]
            # annotations := ((x1+x2)/2,(y1+y2)/2, rank(block.data))
            [x1,x2,x2,x1,x1],[y1,y1,y2,y2,y1]
        end
    end
end
