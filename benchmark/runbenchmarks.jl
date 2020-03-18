using PkgBenchmark
using HierarchicalMatrices

function postprocess(results)
    dir = @__DIR__
    fname = results.commit
    writeresults(joinpath(dir,"results",fname),results)
    export_markdown(joinpath(dir,"results",fname*".md"),results)
end

results   = benchmarkpkg("HierarchicalMatrices")
postprocess(results)
