using BenchmarkTools
using PkgBenchmark
using LinearAlgebra
using HierarchicalMatrices
using HierarchicalMatrices: PartialACA
using HierarchicalMatrices.Clusters: CardinalitySplitter, GeometricMinimalSplitter, AdmissibilityStandard

SUITE                        = BenchmarkGroup()
SUITE["HMatrix"]             = BenchmarkGroup(["hmatrix","hmat"])
SUITE["HMatrix"]["assembly"] = BenchmarkGroup(["assembly","aca"])
SUITE["HMatrix"]["assembly"] = BenchmarkGroup(["assembly","aca"])
SUITE["HMatrix"]["gemv"]     = BenchmarkGroup(["assembly","aca"])

rtol               = 1e-4
N                  = 20000
data               = Geometry.points_on_cylinder(N,1,3/sqrt(N))
f(x,y)::ComplexF64 = x==y ? 0.0 : exp(im*LinearAlgebra.norm(x-y))/LinearAlgebra.norm(x-y)
M                  = LazyMatrix(f,data, data)

splitters   = [Clusters.CardinalitySplitter(),Clusters.GeometricMinimalSplitter()]
admissibles = [Clusters.AdmissibilityStandard(eta) for eta =2:2]
compressors = [HierarchicalMatrices.PartialACA(rtol=rtol)]

for spl in splitters
    clt = Clusters.ClusterTree(data,spl)
    for adm in admissibles
        bclt = Clusters.BlockClusterTree(clt,clt,adm)
        for comp in compressors
            # assembly
            SUITE["HMatrix"]["assembly"]["serial",string(spl),string(adm),string(comp)]  = @benchmarkable $HMatrix($M,$bclt,$comp)
            # SUITE["HMatrix"]["assembly"]["threads",string(spl),string(adm),string(comp)] = @benchmarkable $HMatrix($M,$bclt,$comp;spawn=true)
            # gemv
            x = rand(ComplexF64,N)
            y = similar(x)
            H = HMatrix(M,bclt,comp)
            SUITE["HMatrix"]["gemv"]["serial",string(spl),string(adm),string(comp)] = @benchmarkable mul!($y,$H,$x,1,0)
        end
    end
end
