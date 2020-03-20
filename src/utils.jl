const Maybe{T}         = Union{Tuple{},T}

function debug(flag=true)
    if flag
        @eval ENV["JULIA_DEBUG"] = "HierarchicalMatrices"
    else
        @eval ENV["JULIA_DEBUG"] = ""
    end
end

waitall(t) = for ti in t; wait(ti); end

function LeavesChannel(h::T,args...;kwargs...) where {T}
    chnl = Channel{T}(args...;kwargs...) do chnl
        _collect_leaves(h,chnl)
    end
    return chnl
end

function _collect_leaves(tree,leaves)
    if isleaf(tree)
        put!(leaves,tree)
    else
        for child in getchildren(tree)
            _collect_leaves(child,leaves)
        end
    end
    return leaves
end
