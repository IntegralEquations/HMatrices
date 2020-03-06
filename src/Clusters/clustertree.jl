"""
    ClusterTree{T,S}

Tree structure used to cluster data in of type `T` into containers of type `S`.
"""
mutable struct ClusterTree{T,S}
    data::T
    perm::Vector{Int}
    bounding_box::S
    index_range::UnitRange{Int}
    children::Maybe{Vector{ClusterTree{T,S}}}
    parent::Maybe{ClusterTree{T,S}}
end

AbstractTrees.children(clt::ClusterTree) = getchildren(clt) #interface

getchildren(clt::ClusterTree) = clt.children
getparent(clt::ClusterTree)   = clt.parent
getdata(clt::ClusterTree)     = clt.data
getperm(clt::ClusterTree)     = clt.perm
setchildren!(clt::ClusterTree,children) = (clt.children = children)
setparent!(clt::ClusterTree,parent)     = (clt.parent   = parent)
setdata!(clt::ClusterTree,data)     = (clt.data = data)

isleaf(clt::ClusterTree) = getchildren(clt)  === ()
isroot(clt::ClusterTree) = getparent(clt) === ()

diameter(node::ClusterTree)                         = diameter(node.bounding_box)
distance(node1::ClusterTree,node2::ClusterTree)     = distance(node1.bounding_box, node2.bounding_box)

bounding_box(clt::ClusterTree) = clt.bounding_box

Base.length(node::ClusterTree) = length(node.index_range)
Base.range(node::ClusterTree)  = node.index_range

"""
    ClusterTree(data,splitter)

Construct a `ClusterTree` from the  given `data`.
"""
function ClusterTree(data, splitter=CardinalitySplitter(); reorder=true)
    if reorder
        @info "Input data modified upon construction of ClusterTree"
    else
        data = copy(data)
    end
    n_el    = length(data)
    indices   = collect(1:n_el)
    #build the root, then recurse
    bbox    = bounding_box(data)
    root    = ClusterTree(data,indices,bbox,1:n_el,(),())
    _build_cluster_tree!(root,splitter)
    return root
end

function _build_cluster_tree!(current_node,splitter)
    if should_split(current_node,splitter)
        children          = split!(current_node,splitter)
        setchildren!(current_node,children)
        for child in children
            setparent!(child,current_node)
            _build_cluster_tree!(child,splitter)
        end
    end
end

function Base.show(io::IO,tree::ClusterTree)
    print(io,"ClusterTree with $(length(tree)) elements")
end


################################################################################
## PLOTTING
################################################################################
"""
    plot(tree::ClusterTree,args...)

Plot the point could and the bounding boxes at the leaves of the tree
"""
plot(tree::ClusterTree,args...) = ()

RecipesBase.@recipe function f(tree::ClusterTree)
    legend := false
    grid   --> false
    aspect_ratio --> :equal
    # plot points
    RecipesBase.@series begin
        seriestype := :scatter
        markersize := 2
        tree.data # assumes data Plots knows how  to plot  data
    end
    # plot bounding boxes
    for leaf in AbstractTrees.Leaves(tree)
        RecipesBase.@series begin
            linestyle --> :solid
            color  --> :black
            leaf.bounding_box
        end
    end
end
