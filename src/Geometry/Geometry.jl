module Geometry

import LinearAlgebra
import RecipesBase

import ..Interfaces: container, split, diameter, distance

export Point, HyperRectangle

include("point.jl")
include("utils.jl")
include("hyperrectangle.jl")

end#module
