"""
    SplitStrategy

Splitting strategy used in the creation of a clustertree.

# Options:
- `GeometricSplit` : split each hyper-rectangle in half along the largest axis
- `GeometricMinimalSplit` : like `GeometricSplit`, but shrink to bounding box after each split
- `MedianSplit` : use the median along the largest axis as a splitting point

"""
@enum SplitStrategy  GeometricSplit GeometricMinimalSplit MedianSplit

"""
    ClusterTreeOptions

Specify the options to be used in the  construction of a cluster tree.

# Fields:
- `split_strategy::SplitStrategy`
- `elements_per_leaf::Int`
- `reorder::Bool` : reorder the input points in-place or a make a copy

See also [`SplitStrategy`](@ref)
"""
struct ClusterTreeOptions
    split_strategy::SplitStrategy
    elements_per_leaf::Int
    reorder::Bool
end
ClusterTreeOptions() = ClusterTreeOptions(MedianSplit,100,true)

"""
    ClusterTree{N,T}

Tree structure used to cluster data in the form of `Vector{T}`.
"""
mutable struct ClusterTree{N,T}
    data::Vector{T}
    perm::Vector{Int}
    bounding_box::HyperRectangle{N,Float64}
    index_range::UnitRange{Int}
    children::Maybe{Vector{ClusterTree{N,T}}}
    parent::Maybe{ClusterTree{N,T}}
end

Base.length(node::ClusterTree) = length(node.index_range)
Base.range(node::ClusterTree) = node.index_range
children(clt::ClusterTree) = clt.children
parent(clt::ClusterTree) = clt.parent
isleaf(clt::ClusterTree) = children(clt) === ()
isroot(clt::ClusterTree) = parent(clt) === ()
data(clt::ClusterTree)   = clt.data

function Base.show(io::IO,tree::ClusterTree)
    print(io,"ClusterTree with $(length(tree)) elements")
end

"""
    ClusterTree(data::Vector,options::ClusterTreeOptions)

Construct a `ClusterTree` from the  given `data`.

A typical example is `data::Vector{Point{N,T}}`. The points are then grouped
into bounding boxes of type `HyperRectangle`.

Note that depending on the `reorder` field of ClusterTreeOptions, the passed
data will be modified in-place, therefore mutating the input.

See also: [`ClusterTreeOptions`](@ref)
"""
function ClusterTree(data::Vector, options::ClusterTreeOptions=ClusterTreeOptions())
    if options.reorder
        @info "Input data modified upon construction of ClusterTree"
    else
        data = copy(data)
    end
    n_el    = length(data)
    indices   = collect(1:n_el)
    #build the root, then recurse
    bbox    = Geometry.bounding_box(data)
    root    = ClusterTree(data,indices,bbox,1:n_el,(),())
    _build_cluster_tree!(root,options,indices)
    return root
end

function _build_cluster_tree!(current_node,options,indices)
    if length(current_node) > options.elements_per_leaf
        children          = split(current_node,options,indices)
        current_node.children = children
        for child in children
            child.parent = current_node
            _build_cluster_tree!(child,options,indices)
        end
    end
end

################################################################################
## SPLIT
################################################################################
"""
    split(node::ClusterTree,options)

Divide and redistribute the points of a clustertree node according to the
strategy specified in `options.split_strategy::SplitStrategy`

See also: [`SplitStrategy`](@ref)
"""
function split(node::ClusterTree,options,indices)
    split_strategy = options.split_strategy
    if split_strategy == GeometricSplit
        return _geometric_split(node,indices)
    elseif split_strategy == GeometricMinimalSplit
        return _geometric_minimal_split(node,indices)
    elseif split_strategy == MedianSplit
        return _cardinality_split(node,indices)
    end
end

function _geometric_split(node::ClusterTree,indices)
    rec                 = node.bounding_box
    index_range         = node.index_range
    data                = node.data
    wmax, imax          = findmax(rec.high_corner - rec.low_corner)
    left_rec, right_rec = Geometry.split(rec, imax, rec.low_corner[imax]+wmax/2)
    perm_idxs           = Vector{Int}(undef,length(node))
    npts_left           = 0
    npts_right          = 0
    #sort the points into left and right rectangle
    for i in node.index_range
        pt = data[i]
        if pt in left_rec
            npts_left += 1
            perm_idxs[npts_left] = i
        else
            perm_idxs[length(node)-npts_right] = i
            npts_right += 1
        end
    end
    indices[index_range] = indices[perm_idxs]
    data[index_range]    = data[perm_idxs] # reorders the global index set
    left_index_range     = index_range.start:(index_range.start)+npts_left-1
    right_index_range    = (index_range.start+npts_left):index_range.stop
    return [ClusterTree(data, indices, left_rec,  left_index_range,  (), node),
           ClusterTree(data, indices, right_rec, right_index_range, (), node)]
end

function _geometric_minimal_split(node::ClusterTree,indices)
    rec                 = node.bounding_box
    index_range         = node.index_range
    data                = node.data
    wmax, imax          = findmax(rec.high_corner - rec.low_corner)
    left_rec, right_rec = Geometry.split(rec, imax, rec.low_corner[imax]+wmax/2)
    perm_idxs           = Vector{Int}(undef,length(node))
    npts_left           = 0
    npts_right          = 0
    #sort the points into left and right rectangle
    for i in node.index_range
        pt = data[i]
        if pt in left_rec
            npts_left += 1
            perm_idxs[npts_left] = i
        else
            perm_idxs[length(node)-npts_right] = i
            npts_right += 1
        end
    end
    indices[index_range] = indices[perm_idxs]
    data[index_range]    = data[perm_idxs] # reorders the global index set
    left_index_range     = index_range.start:(index_range.start)+npts_left-1
    right_index_range    = (index_range.start+npts_left):index_range.stop
    # shrink to minimal bounding box before building the next node
    left_rec  = Geometry.bounding_box(data[left_index_range])
    right_rec = Geometry.bounding_box(data[right_index_range])
    return [ClusterTree(data, indices, left_rec,  left_index_range,  (), node),
           ClusterTree(data, indices, right_rec, right_index_range, (), node)]
end

function _cardinality_split(node::ClusterTree,indices)
    rec                 = node.bounding_box
    index_range         = node.index_range
    data                = node.data
    wmax, imax          = findmax(rec.high_corner - rec.low_corner)
    med                 = median(centroid(data[n])[imax] for n in index_range) # the median along largest axis `imax`
    left_rec, right_rec = split(rec, imax, med)
    perm_idxs           = Vector{Int}(undef,length(node))
    npts_left           = 0
    npts_right          = 0
    #sort the points into left and right rectangle
    for i in node.index_range
        pt = centroid(data[i])
        if pt in left_rec
            npts_left += 1
            perm_idxs[npts_left] = i
        else
            perm_idxs[length(node)-npts_right] = i
            npts_right += 1
        end
    end
    indices[index_range] = indices[perm_idxs]
    data[index_range]    = data[perm_idxs] # reorders the global index set
    left_index_range     = index_range.start:(index_range.start)+npts_left-1
    right_index_range    = (index_range.start+npts_left):index_range.stop
    # shrink to minimal bounding box before building the next node
    left_rec  = bounding_box(data[left_index_range])
    right_rec = bounding_box(data[right_index_range])
    return [ClusterTree(data, indices, left_rec,  left_index_range,  (), node),
           ClusterTree(data, indices, right_rec, right_index_range, (), node)]
end

diameter(node::ClusterTree)                         = Geometry.diameter(node.bounding_box)
distance(node1::ClusterTree,node2::ClusterTree)     = distance(node1.bounding_box, node2.bounding_box)

################################################################################
## PLOTTING
################################################################################
"""
    plot(tree::ClusterTree,args...)

Plot the point could and the bounding boxes at the leaves of the tree
"""
plot(tree::ClusterTree,args...) = ()

RecipesBase.@recipe function f(tree::ClusterTree{N}) where {N}
    @assert N<=3 "only 1,2, and 3 dimensional plots available"

    legend := false
    grid   --> false
    aspect_ratio --> :equal
    # plot points
    RecipesBase.@series begin
        seriestype := :scatter
        markersize := 2
        tree.data
    end

    # plot bounding boxes
    for leaf in Leaves(tree)
        RecipesBase.@series begin
            linestyle --> :solid
            color  --> :black
            leaf.bounding_box
        end
    end
end
