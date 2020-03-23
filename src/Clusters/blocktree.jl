"""
    BlockTree{T,S}

Represents a block tree constructed from two `ClusterTree` objects.
"""
mutable struct BlockTree{T,S}
    row_cluster::T
    col_cluster::S
    admissible::Bool
    children::Maybe{Matrix{BlockTree{T,S}}}
    parent::Maybe{BlockTree{T,S}}
end

Base.size(tree::BlockTree)   = (length(tree.row_cluster), length(tree.col_cluster))
Base.size(tree::BlockTree,i) = size(tree)[i]
Base.length(tree::BlockTree) = prod(size(tree))

rowrange(block::BlockTree)   = range(rowcluster(block))
colrange(block::BlockTree)   = range(colcluster(block))

pivot(block::BlockTree) = (rowrange(block).start,colrange(block).start)

AbstractTrees.children(clt::BlockTree) = clt.children #interface

getchildren(bclt::BlockTree)           = bclt.children
getparent(bclt::BlockTree)             = bclt.parent
setchildren!(bclt::BlockTree,children) = (bclt.children = children)
setparent!(bclt::BlockTree,parent)     = (bclt.parent   = parent)

isleaf(clt::BlockTree)       = getchildren(clt) === ()
isroot(clt::BlockTree)       = getparent(clt)   === ()
isadmissible(clt::BlockTree) = clt.admissible

rowcluster(bclt::BlockTree) = bclt.row_cluster
colcluster(bclt::BlockTree) = bclt.col_cluster

"""
    BlockTree(Xtree::ClusterTree,YTree::ClusterTree,admissible_fun=admissible_standard)

Construct a `BlockTree`, and assign to each node a value `admissible` depending on `admissible_fun`

The following signature: `admissible_fun(::BlockTree) --> Bool`
"""
function BlockTree(row_cluster::ClusterTree, col_cluster::ClusterTree, adm_fun=StrongAdmissibilityStd())
    #build root
    root        = BlockTree(row_cluster,col_cluster,false,(),())
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
            block_children     = [BlockTree(r,c,false,(),current_node) for r in row_children, c in col_children]
            setchildren!(current_node,block_children)
            for child in block_children
                _build_block_tree!(adm_fun,child)
            end
        end
    end
end

function Base.show(io::IO,tree::BlockTree)
    print(io,"BlockTree with $(size(tree)) elements")
end

################################################################################
 ## PLOTTING
################################################################################
@recipe function f(tree::BlockTree,filter=(x)->isleaf(x))
    legend --> false
    grid   --> false
    # aspect_ratio --> :equal
    yflip  := true
    seriestype := :shape
    linecolor  --> :black
    # title  := "compression rate = $(compression_rate(hmat))"
    # all leaves
    for block in Iterators.filter(filter,PostOrderDFS(tree))
        @series begin
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
