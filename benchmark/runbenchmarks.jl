using PkgBenchmark
using HierarchicalMatrices
using LinearAlgebra

function postprocess(results,fname=nothing)
    dir = @__DIR__
    fname === nothing && (fname = results.commit)
    writeresults(joinpath(dir,"results",fname),results)
    export_markdown(joinpath(dir,"results",fname*".md"),results)
end

env = Dict("JULIA_NUM_THREADS" =>2,
           "OPEN_BLAS_NUM_THREADS" => 1
           )

config = BenchmarkConfig(env = env,
                         juliacmd = `julia -O3 --check-bounds=no`)

dir       = @__DIR__
fname     = "assembly_benchmark.jl"
# fname     = "gemv_benchmark.jl"
script    = joinpath(dir,fname)
retune    = false
results   = benchmarkpkg("HierarchicalMatrices",config;
                         retune=retune,
                         script=script)
postprocess(results)
