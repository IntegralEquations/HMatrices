abstract type  AbstractHMatrix{T} <: AbstractMatrix{T} end

function Base.getindex(H::AbstractHMatrix,i::Int,j::Int)
    @debug "using `getindex(H::AbstractHMatrix,i::Int,j::Int)`"
    shift = pivot(H) .-1
    _getindex(H,i+shift[1],j+shift[2])
end

function _getindex(H,i,j)
    (i ∈ rowrange(H)) && (j ∈ colrange(H)) || throw(BoundsError(H,(i,j)))
    out = zero(eltype(H))
    if hasdata(H)
        il,jl = idx_global_to_local(i,j,H)
        out  += getdata(H)[il,jl]
    end
    for child in getchildren(H)
        if (i ∈ rowrange(child)) && (j ∈ colrange(child))
            out += _getindex(child,i,j)
        end
    end
    return out
end

Base.size(H::AbstractHMatrix) = length(rowrange(H)), length(colrange(H))

rowrange(H::AbstractHMatrix)         = H.rowrange
colrange(H::AbstractHMatrix)         = H.colrange
pivot(H::AbstractHMatrix)            = (rowrange(H).start,colrange(H).start)
getchildren(H::AbstractHMatrix)      = H.children
getchildren(H::AbstractHMatrix,args...)   = getindex(H.children,args...)
setchildren!(H::AbstractHMatrix,chd) = (H.children = chd)
getparent(H::AbstractHMatrix)        = H.parent
setparent!(H::AbstractHMatrix,par)   = (H.parent   = par)
getdata(H::AbstractHMatrix)          = H.data
setdata!(H::AbstractHMatrix,data)    = (H.data     = data)
isadmissible(H::AbstractHMatrix)     = H.admissible

#indexing by block
struct BlockIndex{T}
    indices::T
end
block(args...) = BlockIndex(args)
blocksize(H::AbstractHMatrix,args...) = size(getchildren(H),args...)
getblock(H::AbstractHMatrix,args...)  = getindex(getchildren(H),args...)
Base.getindex(H::AbstractHMatrix,block::BlockIndex) = getblock(H,block.indices...)

idx_global_to_local(I,J,H::AbstractHMatrix) = (I,J) .- pivot(H) .+ 1
isleaf(H::AbstractHMatrix)                   = getchildren(H) === ()
isroot(H::AbstractHMatrix)                   = getparent(H) === ()
hasdata(H::AbstractHMatrix)                  = getdata(H) !== ()

function Base.zero(H::AbstractHMatrix)
    T = typeof(H)
    if !hasdata(H)
        H0 = T(rowrange(H),colrange(H),isadmissible(H),(),(),())
    else
        H0 = T(rowrange(H),colrange(H),isadmissible(H),zero(getdata(H)),(),())
    end
    if !isleaf(H)
        children = zero.(getchildren(H))
        setchildren!(H0,children)
        map(x->setparent!(x,H0),children)
    end
    return H0
end

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
    for block in PreOrderDFS(H)
        rel_size = length(block)/length(H)
        if hasdata(block)
            c += compression_rate(getdata(block))*rel_size
        end
    end
    return c
end
compression_rate(::Matrix) = 1

# Interface to AbstractTrees
children(H::AbstractHMatrix) = getchildren(H)
Base.eltype(::Type{<:TreeIterator{T}}) where {T<:AbstractHMatrix}         = T
Base.IteratorEltype(::Type{<:TreeIterator{T}}) where {T<:AbstractHMatrix} = Base.HasEltype()

function Base.show(io::IO,hmat::AbstractHMatrix)
    print(io,"hmatrix with range ($(rowrange(hmat))) × ($(colrange(hmat)))")
end

################################################################################
## adjoint
################################################################################
const AdjHMatrix = Adjoint{<:Any,<:AbstractHMatrix}
hasdata(adjH::AdjHMatrix) = hasdata(adjH.parent)
getdata(adjH::AdjHMatrix) = adjoint(getdata(adjH.parent))
getchildren(adjH::AdjHMatrix) = isleaf(adjH.parent) ? () : adjoint(getchildren(adjH.parent)) 
pivot(adjH::AdjHMatrix) = reverse(pivot(adjH.parent))
rowrange(adjH::AdjHMatrix) = colrange(adjH.parent)
colrange(adjH::AdjHMatrix) = rowrange(adjH.parent)
isleaf(adjH::AdjHMatrix) = isleaf(adjH.parent)

function Base.show(io::IO,hmat::AdjHMatrix)
    print(io,"hmatrix with range ($(rowrange(hmat))) × ($(colrange(hmat)))")
end

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
