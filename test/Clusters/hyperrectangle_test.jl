using SafeTestsets

@safetestset "HyperRectangle" begin
    using HierarchicalMatrices.Clusters
    using HierarchicalMatrices.Clusters: HyperRectangle, split, diameter, radius, bounding_box, center
    low_corner  = Point(0.0,0.0)
    high_corner = Point(1.0,2.0)
    mid         = (low_corner + high_corner)/2
    rec = HyperRectangle(low_corner,high_corner)
    @test mid == center(rec)
    @test mid ∈ rec
    @test !in(high_corner + Point(1,1),rec)
    rec1, rec2 = split(rec)
    @test low_corner ∈ rec1
    @test high_corner ∈ rec2
    @test !(low_corner ∈ rec2)
    @test !(high_corner ∈ rec1)
    @test diameter(rec) == sqrt(1^2 + 2^2)
    @test radius(rec) == sqrt(1^2 + 2^2)/2
    # bbox
    pts = Point{2,Float64}[]
    for x=-1:0.1:1
        for y=-1:0.1:1
            push!(pts,Point(x,y))
        end
    end
    @test bounding_box(pts) == HyperRectangle(Point(-1.,-1),Point(1,1.))
end
