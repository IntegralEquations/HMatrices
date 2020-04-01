abstract type  AbstractHMatrix{T} <: AbstractMatrix{T} end

function Base.getindex(H::AbstractHMatrix,i::Int,j::Int)
    (i ∈ rowrange(H)) && (j ∈ colrange(H)) || throw(BoundsError(H,(i,j)))
    out = zero(eltype(H))
    if hasdata(H)
        il,jl = idx_global_to_local(i,j,H)
        out  += getdata(H)[il,jl]
    end
    for child in getchildren(H)
        if (i ∈ rowrange(child)) && (j ∈ colrange(child))
            out += getindex(child,i,j)
        end
    end
    return out
end

Base.size(H::AbstractHMatrix) = length(rowrange(H)), length(colrange(H))

rowrange(H::AbstractHMatrix)         = H.rowrange
colrange(H::AbstractHMatrix)         = H.colrange
pivot(H::AbstractHMatrix)            = (rowrange(H).start,colrange(H).start)
getchildren(H::AbstractHMatrix)      = H.children
setchildren!(H::AbstractHMatrix,chd) = (H.children = chd)
getparent(H::AbstractHMatrix)        = H.parent
setparent!(H::AbstractHMatrix,par)   = (H.parent   = par)
getdata(H::AbstractHMatrix)          = H.data
setdata!(H::AbstractHMatrix,data)    = (H.data     = data)
isadmissible(H::AbstractHMatrix)     = H.admissible

idx_global_to_local(I,J,H::AbstractHMatrix) = (I,J) .- pivot(H) .+ 1
isleaf(H::AbstractHMatrix)                   = getchildren(H) === ()
isroot(H::AbstractHMatrix)                   = getparent(H) === ()
hasdata(H::AbstractHMatrix)                  = getdata(H) !== ()

function Matrix(H::AbstractHMatrix)
    M = zeros(eltype(H),size(H)...)
    shift = pivot(H) .- 1
    for block in PreOrderDFS(H)
        hasdata(block) || continue
        irange = rowrange(block) .- shift[1]
        jrange = colrange(block) .- shift[2]
        M[irange,jrange] += Matrix(block.data)
    end
    return M
end

function compression_rate(H::AbstractHMatrix)
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
children(H::AbstractHMatrix) = getchildren(H)
Base.eltype(::Type{<:TreeIterator{T}}) where {T<:AbstractHMatrix}         = T
Base.IteratorEltype(::Type{<:TreeIterator{T}}) where {T<:AbstractHMatrix} = Base.HasEltype()

################################################################################
## Plot recipes
################################################################################
@recipe function f(hmat::AbstractHMatrix)
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
