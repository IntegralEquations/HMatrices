# a place to declare the functions gluing the modules together

module Interfaces

# Geometry <--> Clusters
function bounding_box end
function split end

# trees
function getchildren end
function getparent end
function setchildren end
function setparent end

end#module
