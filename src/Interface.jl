# a place to declare the functions gluing the modules together

module Interface

# trees
function getchildren end
function getparent end
function setchildren end
function setparent end
function isleaf end
function isroot end

# Cluster <--> HMatrices
function rowrange end
function colrange end
function isadmissible end

end#module
