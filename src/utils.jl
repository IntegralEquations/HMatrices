function debug(flag=true)
    if flag
        @eval ENV["JULIA_DEBUG"] = "HMatrices"
    else
        @eval ENV["JULIA_DEBUG"] = ""
    end
end

waitall(t) = for ti in t; wait(ti); end

function depth(node,acc=0)
    if isroot(node)
        return acc
    else
        depth(getparent(node),acc+1)
    end
end

# function points_on_cylinder(n, radius, step)
#     result          = Vector{Point{3,Float64}}(undef,n)
#     length          = 2*π*radius
#     pointsPerCircle = length/step
#     angleStep       = 2*π/pointsPerCircle
#     for i=0:n-1
#         x = radius * cos(angleStep*i)
#         y = radius * sin(angleStep*i)
#         z = step*i/pointsPerCircle
#         result[i+1] = Point(x,y,z)
#     end
#     return result
# end
