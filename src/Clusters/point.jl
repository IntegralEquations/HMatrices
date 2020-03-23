using GeometryTypes: Point

#TODO: make native point type

# Base.min(a::Point,b::Point) = Point(min.(a,b))
# Base.max(a::Point,b::Point) = Point(max.(a,b))

# struct Point{N,T}
#     coords::NTuple{N,T}
# end

# Point(args...) = Point(promote(args...))

# Base.promote_rule(::Type{Point{N,T}},::Type{Point{N,S}}) where {N,T,S} = Point{N,promote_type(T,S)}

# Base.length(pt::Point{N}) where {N} = N
# Base.iterate(pt::Point,args...) = iterate(pt.coords,args...)

# Base.convert(::Type{T},x::NTuple) where {T <: Point} = T(x)
# Base.convert(::Type{T},x::Point) where {T <: Point} = T(x.coords)

# for op in (:+, :-, :min, :max)
#     @eval Base.$op(pt1::Point{N},pt2::Point{N}) where {N} = Point($op.(pt1.coords,pt2.coords))
# end

# Base.getindex(pt::Point,i) = getindex(pt.coords,i)
