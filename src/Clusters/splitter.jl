abstract type AbstractSplitter end

@Base.kwdef struct GeometricSplitter <: AbstractSplitter
    nmax::Int=Parameters.nmax
end

@Base.kwdef struct GeometricMinimalSplitter <: AbstractSplitter
    nmax::Int=Parameters.nmax
end

@Base.kwdef struct CardinalitySplitter <: AbstractSplitter
    nmax::Int=Parameters.nmax
end

@Base.kwdef struct NestedDissectionSplitter <: AbstractSplitter
    nmax::Int=Parameters.nmax
end

should_split(node::ClusterTree,splitter::AbstractSplitter) = length(node) > splitter.nmax

function _binary_split(cluster::ClusterTree,dir,pos,shrink=true)
    perm                = cluster.perm
    rec                 = container(cluster)
    index_range         = cluster.index_range
    data                = getdata(cluster)
    left_rec, right_rec = split(rec, dir, pos)
    perm_idxs           = Vector{Int}(undef,length(cluster))
    npts_left           = 0
    npts_right          = 0
    #sort the points into left and right rectangle
    for i in cluster.index_range
        pt = data[i]
        if pt in left_rec
            npts_left += 1
            perm_idxs[npts_left] = i
        else
            perm_idxs[length(cluster)-npts_right] = i
            npts_right += 1
        end
    end
    perm[index_range]     = perm[perm_idxs]
    data[index_range]     = data[perm_idxs] # reorders the global index set
    left_index_range      = index_range.start:(index_range.start)+npts_left-1
    right_index_range     = (index_range.start+npts_left):index_range.stop
    if shrink
        left_rec   = container(data[left_index_range])
        right_rec  = container(data[right_index_range])
    end
    clt1 = ClusterTree(data, perm, left_rec,  left_index_range,  (), cluster)
    clt2 = ClusterTree(data, perm, right_rec, right_index_range, (), cluster)
    return [clt1, clt2]
end

function split!(cluster::ClusterTree,splitter::GeometricSplitter)
    rec  = container(cluster)
    wmax, imax          = findmax(rec.high_corner - rec.low_corner)
    return _binary_split(cluster, imax, rec.low_corner[imax]+wmax/2, false)
end

function split!(cluster::ClusterTree,splitter::GeometricMinimalSplitter)
    rec  = container(cluster)
    wmax, imax          = findmax(rec.high_corner - rec.low_corner)
    return _binary_split(cluster, imax, rec.low_corner[imax]+wmax/2)
end

function split!(cluster::ClusterTree,splitter::CardinalitySplitter)
    data                = getdata(cluster)
    rec                 = container(cluster)
    index_range         = cluster.index_range
    wmax, imax          = findmax(rec.high_corner - rec.low_corner)
    med                 = median(data[n][imax] for n in index_range) # the median along largest axis `imax`
    return _binary_split(cluster, imax, med)
end

function split!(cluster::ClusterTree,splitter::NestedDissectionSplitter)
    #TODO
    # clt1, clt2 = split(cluster,GeometricMinimalSplitter(splitter.nmax))
    # S          = getconnectivity(getdata(cluster))
    # perm       = cluster.perm
    # tmp        = Int[]
    # for i in clt1.index_range
    #     for j in clt2.index_range
    #         if S[perm[i],perm[j]] != 0
    #             push!(tmp,j)
    #             push!(tmp,i)
    #         end
    #     end
    # end
    # # sep = Cluster(data,perm,bbox,sep_index_range,(),cluster)
    # return [clt1, clt2, tmp]
end
