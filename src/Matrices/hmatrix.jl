mutable struct HMatrix{S,F,T} <: AbstractHierarchicalMatrix{T}
    rowrange::UnitRange{Int}
    colrange::UnitRange{Int}
    data::Maybe{Union{S,F}}
    children::Maybe{Matrix{HMatrix{S,F,T}}}
    parent::Maybe{HMatrix{S,F,T}}
    function HMatrix{S,F,T}(irange,jrange,data,children,parent) where {S,F,T}
        if (data !== ())
            @assert (length(irange),length(jrange)) === size(data) "$(length(irange)),$(length(jrange)) != $(size(data))"
        end
        new{S,F,T}(irange,jrange,data,children,parent)
    end
end

HMatrix{S,F,T}(block::BlockClusterTree) where  {S,F,T} = HMatrix{S,F,T}(rowrange(block),colrange(block),(),(),())

function HMatrix{S,F,T}(K,block,comp) where {S,F,T}
    #TODO: add prehook?
    # initialize an empty root
    hmat = HMatrix{S,F,T}(block)
    # recursion
    _build_hmatrix!(hmat,block,K,comp) # recursion function
    #TODO: add posthook to e.g. coarsen or recompress immediately  after construction
    return hmat
end

function _build_hmatrix!(hmat,block,K,comp)
    T = eltype(hmat)
    S = sparsetype(hmat)
    F = densetype(hmat)
    if !isleaf(block)
        children = HMatrix{S,F,T}.(getchildren(block)) # initalize empty sub-blocks
        setchildren!(hmat,children)
        map(x->setparent!(x,hmat),children)
        for (hmat_child,block_child) in zip(getchildren(hmat),getchildren(block))
            _build_hmatrix!(hmat_child,block_child,K,comp)
        end
    elseif !isadmissible(block)
        data = [K[i,j] for i in rowrange(block), j in colrange(block)]
        setdata!(hmat,data)
    else
        data = compress(K,block,comp)
        setdata!(hmat,data)
    end
    return hmat
end

function HMatrix(K,block,comp)
    T = eltype(K)
    S = Base.promote_op(compress,typeof(K),typeof(block),typeof(comp))
    F = Matrix{T}
    HMatrix{S,F,T}(K,block,comp)
end

sparsetype(h::HMatrix{S}) where {S} = S
densetype(h::HMatrix{_,F}) where {_,F} = F


function Base.show(io::IO,hmat::HMatrix)
    print(io,"hmatrix with range ($(rowrange(hmat))) × ($(colrange(hmat)))")
end
