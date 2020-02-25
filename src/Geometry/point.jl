struct Point{N,T}
    coords::NTuple{N,T}
end
Base.size(pt::Point{N,T}) where {N,T} = (N,)

Point(args...) = Point(promote(args...))

for op in (:+, :-, :min, :max)
    @eval Base.$op(pt1::T,pt2::T) where {T<:Point} = T($op.(pt1.coords,pt2.coords))
end

Base.getindex(pt::Point,i) = getindex(pt.coords,i)
