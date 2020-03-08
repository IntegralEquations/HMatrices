function points_on_cylinder(n, radius, step)
    result          = Vector{Point{3,Float64}}(undef,n)
    length          = 2*π*radius
    pointsPerCircle = length/step
    angleStep       = 2*π/pointsPerCircle
    for i=0:n-1
        x = radius * cos(angleStep*i)
        y = radius * sin(angleStep*i)
        z = step*i/pointsPerCircle
        result[i+1] = Point(x,y,z)
    end
    return result
end
