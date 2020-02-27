"""
    BlockClusterTree{T,S}

Represents a block tree constructed from two `ClusterTree` objects.
"""
mutable struct BlockClusterTree{T,S}
    row_cluster::T
    col_cluster::S
    admissible::Bool
    children::Maybe{Matrix{BlockClusterTree{T,S}}}
    parent::Maybe{BlockClusterTree{T,S}}
end

Base.size(tree::BlockClusterTree)   = (length(tree.row_cluster), length(tree.col_cluster))
Base.size(tree::BlockClusterTree,i) = size(tree)[i]
Base.length(tree::BlockClusterTree) = prod(size(tree))

rowrange(block::BlockClusterTree)   = range(rowcluster(block))
colrange(block::BlockClusterTree)   = range(colcluster(block))

pivot(block::BlockClusterTree) = (rowrange(block).start,colrange(block).start)

AbstractTrees.children(clt::BlockClusterTree) = clt.children #interface

getchildren(bclt::BlockClusterTree)           = bclt.children
getparent(bclt::BlockClusterTree)             = bclt.parent
setchildren!(bclt::BlockClusterTree,children) = (bclt.children = children)
setparent!(bclt::BlockClusterTree,parent)     = (bclt.parent   = parent)

isleaf(clt::BlockClusterTree)       = children(clt) === ()
isroot(clt::BlockClusterTree)       = parent(clt)   === ()
isadmissible(clt::BlockClusterTree) = clt.admissible

rowcluster(bclt::BlockClusterTree) = bclt.row_cluster
colcluster(bclt::BlockClusterTree) = bclt.col_cluster

"""
    BlockClusterTree(Xtree::ClusterTree,YTree::ClusterTree,admissible_fun=admissible_standard)

Construct a `BlockClusterTree`, and assign to each node a value `admissible` depending on `admissible_fun`

The following signature: `admissible_fun(::BlockClusterTree) --> Bool`
"""
function BlockClusterTree(row_cluster::ClusterTree, col_cluster::ClusterTree, adm_fun=AdmissibiltyStandard())
    #build root
    root        = BlockClusterTree(row_cluster,col_cluster,false,(),())
    # recurse
    _build_block_tree!(adm_fun,root)
    return root
end

function _build_block_tree!(adm_fun, current_node::T) where {T}
    if adm_fun(current_node)
        current_node.admissible = true
    else
        current_node.admissible = false
        if !(isleaf(rowcluster(current_node)) || isleaf(colcluster(current_node)))
            row_children       = getchildren(rowcluster(current_node))
            col_children       = getchildren(colcluster(current_node))
            block_children     = [BlockClusterTree(r,c,false,(),current_node) for r in row_children, c in col_children]
            setchildren!(current_node,block_children)
            for child in block_children
                _build_block_tree!(adm_fun,child)
            end
        end
    end
end

function Base.show(io::IO,tree::BlockClusterTree)
    print(io,"BlockClusterTree with $(size(tree)) elements")
end

################################################################################
 ## PLOTTING
################################################################################
RecipesBase.@recipe function f(bclt::BlockClusterTree)
    legend --> false
    grid   --> false
    # aspect_ratio --> :equal
    yflip  := true
    seriestype := :shape
    linecolor  --> :black
    # title  := "compression rate = $(compression_rate(hmat))"
    # all leaves
    for block in AbstractTrees.Leaves(bclt)
        RecipesBase.@series begin
            if isadmissible(block)
                fillcolor  := :blue
                opacity    := 0.5
            else
                fillcolor  := :red
                opacity     = 0.3
            end
            pt1 = pivot(block)
            pt2 = pt1 .+ size(block) .- 1
            y1, y2 = pt1[1],pt2[1]
            x1, x2 = pt1[2],pt2[2]
            [x1,x2,x2,x1,x1],[y1,y1,y2,y2,y1]
        end
    end
end
