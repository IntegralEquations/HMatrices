using HierarchicalMatrices: FlexMatrix
using HierarchicalMatrices: compress, ACA, PartialACA
using BenchmarkTools
using HierarchicalMatrices
using Random

using HierarchicalMatrices.Geometry
using HierarchicalMatrices.Clusters
using LinearAlgebra

Random.seed!(0)

const rtol = 1e-4

N           = 10000
data        = Geometry.points_on_cylinder(N,1,3/sqrt(N))
f(x,y)::ComplexF64 = x==y ? 0.0 : exp(im*LinearAlgebra.norm(x-y))/LinearAlgebra.norm(x-y)
M        = LazyMatrix(f,data, data)

SUITE["assembly"] = BenchmarkGroup()

splitters   = [Clusters.CardinalitySplitter(),Clusters.GeometricMinimalSplitter()]
admissibles = [Clusters.AdmissibilityStandard(eta) for eta =1:1]
compressors = [HierarchicalMatrices.PartialACA(rtol=rtol), HierarchicalMatrices.ACA(rtol=rtol)]
for spl in splitters, adm in admissibles, comp in compressors
    clt                        = ClusterTree(data,spl)
    bclt                       = BlockClusterTree(clt,clt,adm)
    SUITE["assembly"][string(spl),string(adm),string(comp)]  = @benchmarkable $HMatrix($M,$bclt,$comp)
end

