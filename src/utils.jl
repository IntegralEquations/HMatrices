const Maybe{T}         = Union{Tuple{},T}

function debug(flag=true)
    if flag
        @eval ENV["JULIA_DEBUG"] = "HierarchicalMatrices"
    else
        @eval ENV["JULIA_DEBUG"] = ""
    end
end
