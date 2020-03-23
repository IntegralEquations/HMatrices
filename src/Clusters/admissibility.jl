abstract type AbstractAdmissibilityCondition end

(adm::AbstractAdmissibilityCondition)(bclt::BlockTree) = adm(rowcluster(bclt),colcluster(bclt))

Base.@kwdef struct StrongAdmissibilityStd <: AbstractAdmissibilityCondition
    eta::Float64=Parameters.eta
end

function (adm::StrongAdmissibilityStd)(left_node::ClusterTree, right_node::ClusterTree)
    diam_min = minimum(diameter,(left_node,right_node))
    dist     = distance(left_node,right_node)
    return diam_min < adm.eta*dist
end

struct WeakAdmissibilityStd <: AbstractAdmissibilityCondition end
(adm::WeakAdmissibilityStd)(left_node::ClusterTree, right_node::ClusterTree) = distance(left_node,right_node) > 0
