using BenchmarkTools
using PkgBenchmark

SUITE = BenchmarkGroup()

include("assembly_benchmark.jl")
