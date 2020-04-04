abstract type  HierarchicalMatrix{T} <: AbstractMatrix{T} end

mutable struct HMatrix{S,F,T} <: HierarchicalMatrix{T}
    rowrange::UnitRange{Int}
    colrange::UnitRange{Int}
    admissible::Bool
    data::Maybe{Union{S,F}}
    children::Maybe{Matrix{HMatrix{S,F,T}}}
    parent::Maybe{HMatrix{S,F,T}}
    function HMatrix{S,F,T}(irange,jrange,adm,data,children,parent) where {S,F,T}
        if (data !== ())
            @assert (length(irange),length(jrange)) === size(data) "$(length(irange)),$(length(jrange)) != $(size(data))"
        end
        new{S,F,T}(irange,jrange,adm,data,children,parent)
    end
end

HMatrix(args...;kwargs...) = HMatrix(CPU1(),args...;kwargs...)

function Base.getindex(H::HierarchicalMatrix,i::Int,j::Int)
    @debug "using `getindex(H::HierarchicalMatrix,i::Int,j::Int)`"
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

Base.size(H::HierarchicalMatrix) = length(rowrange(H)), length(colrange(H))

rowrange(H::HierarchicalMatrix)         = H.rowrange
colrange(H::HierarchicalMatrix)         = H.colrange
pivot(H::HierarchicalMatrix)            = (rowrange(H).start,colrange(H).start)
getchildren(H::HierarchicalMatrix)      = H.children
getchildren(H::HierarchicalMatrix,args...)   = getindex(H.children,args...)
setchildren!(H::HierarchicalMatrix,chd) = (H.children = chd)
getparent(H::HierarchicalMatrix)        = H.parent
setparent!(H::HierarchicalMatrix,par)   = (H.parent   = par)
getdata(H::HierarchicalMatrix)          = H.data
setdata!(H::HierarchicalMatrix,data)    = (H.data     = data)
isadmissible(H::HierarchicalMatrix)     = H.admissible

#indexing by block
struct BlockIndex{T}
    indices::T
end
block(args...) = BlockIndex(args)
blocksize(H::HierarchicalMatrix,args...) = size(getchildren(H),args...)
getblock(H::HierarchicalMatrix,args...)  = getindex(getchildren(H),args...)
Base.getindex(H::HierarchicalMatrix,block::BlockIndex) = getblock(H,block.indices...)

idx_global_to_local(I,J,H::HierarchicalMatrix) = (I,J) .- pivot(H) .+ 1
isleaf(H::HierarchicalMatrix)                   = getchildren(H) === ()
isroot(H::HierarchicalMatrix)                   = getparent(H) === ()
hasdata(H::HierarchicalMatrix)                  = getdata(H) !== ()

function Base.zero(H::HierarchicalMatrix)
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

function Matrix(H::HierarchicalMatrix)
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

function compression_rate(H::HierarchicalMatrix)
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
children(H::HierarchicalMatrix) = getchildren(H)
Base.eltype(::Type{<:TreeIterator{T}}) where {T<:HierarchicalMatrix}         = T
Base.IteratorEltype(::Type{<:TreeIterator{T}}) where {T<:HierarchicalMatrix} = Base.HasEltype()

function Base.show(io::IO,hmat::HierarchicalMatrix)
    print(io,"hmatrix with range ($(rowrange(hmat))) × ($(colrange(hmat)))")
end

"""
    HMatrix([resource=CPU1()],K,blocktree,comp)

Main constructor for hierarchical matrix, where `K` represents the matrix to be
approximated, `blocktree` encondes the tree structure, and `comp` is a
function/functor which can compress blocks.

It is assumed that `K` supports `getindex(K,i,j)`, that `blocktree` has methods
`getchildren(blocktree)`, `isadmissible(blocktree)`,
`rowrange(blocktree)->UnitRange`, and `colrange(blocktree)->UnitRange` , and
that comp can be called as `comp(K,irange::UnitRange,jrange::UnitRange)` to
produce a compressed version of `K[irange,jrange]`.

An optional first argument can be passed to how to perform the computation (e.g. CPU1, CPUThreads).
"""
function HMatrix(resource::AbstractResource,K,blocktree,comp)
    T = eltype(K)
    S = Base.promote_op(comp,typeof(K),typeof(rowrange(blocktree)),typeof(colrange(blocktree)))
    F = Matrix{T}
    HMatrix{S,F,T}(resource,K,blocktree,comp)
end

function HMatrix{S,F,T}(resource::AbstractResource,K,block,comp) where {S,F,T}
    #TODO: add prehook?
    # initialize the HMatrix structure from block
    hmat = HMatrix{S,F,T}(block)
    # recursion
    assemble!(resource,hmat,K,comp) # recursive function
    #TODO: add posthook to e.g. coarsen or recompress immediately  after construction
    return hmat
end

function HMatrix{S,F,T}(block::BlockTree) where {S,F,T}
    hmat = HMatrix{S,F,T}(rowrange(block),colrange(block),isadmissible(block),(),(),())
    children = HMatrix{S,F,T}.(getchildren(block))
    setchildren!(hmat,children)
    map(x->setparent!(x,hmat),children)
    return hmat
end

function assemble!(resource::CPU1,hmat,K,comp)
    if isleaf(hmat)
        if isadmissible(hmat)
            data = comp(K,rowrange(hmat),colrange(hmat))
            setdata!(hmat,data)
        else
            data = [K[i,j] for i in rowrange(hmat), j in colrange(hmat)]
            setdata!(hmat,data)
        end
    else
        for child in getchildren(hmat)
            assemble!(resource,child,K,comp)
        end
    end
    return hmat
end

################################################################################
## Threaded assemble of hmat
################################################################################
# TODO: implement a filter to control the granularity of the threaded assemble
Base.@kwdef struct ThreadedAssembleOpts{T}
    filter::T= (x)->isleaf(x)
end

function assemble!(resources::CPUThreads,hmat,K,comp)
    settings = resources.settings
    if settings === nothing
        @sync for node in Leaves(hmat)
            @spawn assemble!(CPU1(),node,K,comp)
        end
    else
        @sync for node in Iterators.filter(settings.filter,PostOrderDFS(hmat))
            @spawn assemble!(CPU1(),node,K,comp)
        end
    end
    return hmat
end

sparsetype(h::HMatrix{S}) where {S} = S
densetype(h::HMatrix{_,F}) where {_,F} = F

################################################################################
## Plot recipes
################################################################################
@recipe function f(hmat::HierarchicalMatrix)
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
