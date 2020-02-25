"""
    HyperRectangle{N,T}
Structure representing a hyperrectangle in `N` dimension with datatype `T`,
specified by a `low_corner` and a `high_corner`
"""
struct HyperRectangle{N,T}
    low_corner::Point{N,T}
    high_corner::Point{N,T}
end

Base.:(==)(h1::HyperRectangle, h2::HyperRectangle) = (h1.low_corner == h2.low_corner) && (h1.high_corner == h2.high_corner)
Base.in(point,h::HyperRectangle) = all(h.high_corner .>= point .>= h.low_corner)
Base.eltype(h::HyperRectangle{N,T}) where {N,T} = T
dimension(h::HyperRectangle{N}) where {N} = N

################################################################################
## CONVENIENCE FUNCTIONS FOR HYPERRECTANGLES
################################################################################

"""
    split(rec::HyperRectangle,axis::Int,place)
Split a hyperrectangle in two along the `axis` direction at the  position `place`. Returns the two resulting hyper-rectangles.
"""
function split(rec::HyperRectangle,axis,place)
    N            = dimension(rec)
    T            = eltype(rec)
    high_corner1 = ntuple(n-> n==axis ? place : rec.high_corner[n], N)
    low_corner2  = ntuple(n-> n==axis ? place : rec.low_corner[n], N)
    rec1         = HyperRectangle{N,T}(rec.low_corner, high_corner1)
    rec2         = HyperRectangle{N,T}(low_corner2,rec.high_corner)
    return [rec1, rec2]
end

"""
    split(rec::HyperRectangle,axis)
When no `place` is given, defaults to splitting in the middle of the axis.
"""
function split(rec::HyperRectangle,axis)
    place              = (rec.high_corner[axis] + rec.low_corner[axis])/2
    split(rec,axis,place)
end


"""
    split(rec::HyperRectangle)
When no axis and no place is given, defaults to splitting along the largest axis
"""
function split(rec::HyperRectangle)
    axis = argmax(rec.high_corner .- rec.low_corner)
    split(rec,axis)
end

diameter(cub::HyperRectangle) = LinearAlgebra.norm(cub.high_corner - cub.low_corner,2)

function distance(rec1::HyperRectangle{N},rec2::HyperRectangle{N}) where {N}
    d2 = 0
    for i=1:N
        d2 += max(0,rec1.low_corner[i] - rec2.high_corner[i])^2 + max(0,rec2.low_corner[i] - rec1.high_corner[i])^2
    end
    return sqrt(d2)
end

centroid(x) = x

function bounding_box(data::Vector{<:Point})
    pt_min = minimum(data)
    pt_max = maximum(data)
    return HyperRectangle(pt_min,pt_max)
end
bounding_box(data::Vector) = bounding_box(centroid.(data))

center(rec::HyperRectangle) = (rec.low_corner + rec.high_corner) ./ 2
radius(rec::HyperRectangle) = norm(rec.high_corner-rec.low_corner) ./ 2


################################################################################
## PLOTTING RECIPES
################################################################################
@recipe function f(rec::HyperRectangle{N}) where {N}
    if N == 2
        seriestype := :path
        pt1 = rec.low_corner
        pt2 = rec.high_corner
        x1, x2 = pt1[1],pt2[1]
        y1, y2 = pt1[2],pt2[2]
        [x1,x2,x2,x1,x1],[y1,y1,y2,y2,y1]
    elseif N == 3
        seriestype := :path
        pt1 = rec.low_corner
        pt2 = rec.high_corner
        x1, x2 = pt1[1],pt2[1]
        y1, y2 = pt1[2],pt2[2]
        z1, z2 = pt1[3],pt2[3]
        @series [x1,x2,x2,x1,x1],[y1,y1,y2,y2,y1],[z1,z1,z1,z1,z1]
        @series [x1,x2,x2,x1,x1],[y1,y1,y2,y2,y1],[z2,z2,z2,z2,z2]
        @series [x1,x1],[y1,y1],[z1,z2]
        @series [x2,x2],[y1,y1],[z1,z2]
        @series [x2,x2],[y2,y2],[z1,z2]
        @series [x1,x1],[y2,y2],[z1,z2]
    end
end
