# a place to declare the functions gluing the modules together

module Interfaces

# Geometry <--> Clusters
function container end
function split end
function diameter end
function distance end

# trees
function getchildren end
function getparent end
function setchildren end
function setparent end

end#module
