using BenchmarkTools
using PkgBenchmark
using LinearAlgebra
using HierarchicalMatrices
using HierarchicalMatrices: PartialACA
# using HierarchicalMatrices.Clusters: CardinalitySplitter, GeometricMinimalSplitter, AdmissibilityStandard
using ComputationalResources

const SUITE = BenchmarkGroup()
SUITE["HMatrix"] = BenchmarkGroup(["hmatrix", "hmat"])
SUITE["HMatrix"]["assembly"] = BenchmarkGroup(["assembly", "aca"])

rtol = 1e-4
N = 20000
data = Geometry.points_on_cylinder(N, 1, 3 / sqrt(N))
f(x, y)::ComplexF64 = x == y ? 0.0 :
    exp(im * LinearAlgebra.norm(x - y)) / LinearAlgebra.norm(x - y)
M = LazyMatrix(f, data, data)

spl = Clusters.GeometricMinimalSplitter()
adm = Clusters.AdmissibilityStandard(3)
comp = HierarchicalMatrices.PartialACA(rtol = rtol)

clt = Clusters.ClusterTree(data, spl)
bclt = Clusters.BlockTree(clt, clt, adm)

H = HMatrix(CPUThreads(), M, bclt, comp)

SUITE["HMatrix"]["assembly"]["CPU1"]       = @benchmarkable $HMatrix($M, $bclt, $comp)
SUITE["HMatrix"]["assembly"]["CPUThreads"] = @benchmarkable $HMatrix($(CPUThreads()), $M, $bclt, $comp)
