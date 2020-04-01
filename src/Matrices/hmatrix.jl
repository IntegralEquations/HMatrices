mutable struct HMatrix{S,F,T} <: AbstractHMatrix{T}
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
function HMatrix(resource::AbstractCPU,K,blocktree,comp)
    T = eltype(K)
    S = Base.promote_op(comp,typeof(K),typeof(rowrange(blocktree)),typeof(colrange(blocktree)))
    F = Matrix{T}
    HMatrix{S,F,T}(resource,K,blocktree,comp)
end

function HMatrix{S,F,T}(resource::AbstractCPU,K,block,comp) where {S,F,T}
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

function RkMatrix(H::HMatrix)
    children = getchildren(H)
    if isleaf(H)
        hasdata(H) && (data = getdata(H))
        if data isa RkMatrix
            return data
        else
            return RkMatrix(data)
        end
    else
        B = [RkMatrix(child) for child in getchildren(H)]
        tmp1 = hcat(B[1,1],B[1,2])
        tmp2 = hcat(B[2,1],B[2,2])
        tmp3 = vcat(tmp1,tmp2)
        return tmp3
    end
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

function Base.show(io::IO,hmat::HMatrix)
    print(io,"hmatrix with range ($(rowrange(hmat))) × ($(colrange(hmat)))")
end
