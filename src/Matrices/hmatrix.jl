mutable struct HMatrix{S,F,T} <: AbstractHierarchicalMatrix{T}
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

function HMatrix(K,block,comp)
    T = eltype(K)
    S = Base.promote_op(comp,typeof(K),typeof(rowrange(block)),typeof(colrange(block)))
    F = Matrix{T}
    HMatrix{S,F,T}(K,block,comp)
end

function HMatrix{S,F,T}(K,block,comp) where {S,F,T}
    #TODO: add prehook?
    # initialize the HMatrix structure from block
    hmat = HMatrix{S,F,T}(block)
    # recursion
    assemble!(hmat,K,comp) # recursive function
    #TODO: add posthook to e.g. coarsen or recompress immediately  after construction
    return hmat
end

function assemble!(hmat,K,comp)
    if !isleaf(hmat)
        for child in getchildren(hmat)
            assemble!(child,K,comp)
        end
    elseif !isadmissible(hmat)
        data = [K[i,j] for i in rowrange(hmat), j in colrange(hmat)]
        setdata!(hmat,data)
    else
        data = comp(K,rowrange(hmat),colrange(hmat))
        setdata!(hmat,data)
    end
    return hmat
end

function HMatrix{S,F,T}(block::BlockClusterTree) where {S,F,T}
    hmat = HMatrix{S,F,T}(rowrange(block),colrange(block),isadmissible(block),(),(),())
    children = HMatrix{S,F,T}.(getchildren(block))
    setchildren!(hmat,children)
    map(x->setparent!(x,hmat),children)
    return hmat
end

sparsetype(h::HMatrix{S}) where {S} = S
densetype(h::HMatrix{_,F}) where {_,F} = F

function Base.show(io::IO,hmat::HMatrix)
    print(io,"hmatrix with range ($(rowrange(hmat))) × ($(colrange(hmat)))")
end
