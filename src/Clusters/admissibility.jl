abstract type AbstractAdmissibilityCondition end

(adm::AbstractAdmissibilityCondition)(bclt::BlockClusterTree) = adm(rowcluster(bclt),colcluster(bclt))

Base.@kwdef struct AdmissibilityStandard <: AbstractAdmissibilityCondition
    eta::Float64=2
end

function (adm::AdmissibilityStandard)(left_node::ClusterTree, right_node::ClusterTree)
    diam_min = minimum(diameter,(left_node,right_node))
    dist     = distance(left_node,right_node)
    return diam_min < adm.eta*dist
end
