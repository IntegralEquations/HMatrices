const AdjHMatrix = Adjoint{<:Any,<:HierarchicalMatrix}
hasdata(adjH::AdjHMatrix) = hasdata(adjH.parent)
getdata(adjH::AdjHMatrix) = adjoint(getdata(adjH.parent))
getchildren(adjH::AdjHMatrix) = isleaf(adjH.parent) ? () : adjoint(getchildren(adjH.parent)) 
pivot(adjH::AdjHMatrix) = reverse(pivot(adjH.parent))
rowrange(adjH::AdjHMatrix) = colrange(adjH.parent)
colrange(adjH::AdjHMatrix) = rowrange(adjH.parent)
isleaf(adjH::AdjHMatrix) = isleaf(adjH.parent)

function Base.show(io::IO,hmat::AdjHMatrix)
    print(io,"hmatrix with range ($(rowrange(hmat))) × ($(colrange(hmat)))")
end
