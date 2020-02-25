using SafeTestsets

@safetestset "Point" begin
    using HierarchicalMatrices.Geometry
    pt1 = Point(2,2)
    pt2 = Point(1,3)
    @test pt1 + pt2 == Point(3,5)
    @test pt1 - pt2 == Point(1,-1)
    @test min(pt1,pt2) == Point(1,2)
    @test max(pt1,pt2) == Point(2,3)
    pt1 = Point(2,2,3.)
    pt2 = Point(1,3,4.)
    @test pt1 + pt2    == Point(3,5,7.)
    @test pt1 - pt2    == Point(1,-1,-1.)
    @test min(pt1,pt2) == Point(1,2,3.)
    @test max(pt1,pt2) == Point(2,3,4.)
end

