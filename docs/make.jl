push!(LOAD_PATH,joinpath(@__DIR__, ".."))
using Documenter, HierarchicalMatrices

makedocs(
    modules = [HierarchicalMatrices],
    format = Documenter.HTML(; prettyurls = get(ENV, "CI", nothing) == "true"),
    authors = "Luiz M. Faria",
    sitename = "HierarchicalMatrices.jl",
    pages = Any["index.md"]
    # strict = true,
    # clean = true,
    # checkdocs = :exports,
)

deploydocs(
    repo = "github.com/maltezfaria/HierarchicalMatrices.jl.git",
)
