module Geometry

import LinearAlgebra
import RecipesBase

import ..Interfaces: bounding_box, split

export Point, HyperRectangle

include("point.jl")
include("hyperrectangle.jl")

end#module
