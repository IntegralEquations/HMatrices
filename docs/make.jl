push!(LOAD_PATH,joinpath(@__DIR__, ".."))
using Documenter, HierarchicalMatrices

makedocs(
    modules = [HierarchicalMatrices],
    format = Documenter.HTML(; prettyurls = get(ENV, "CI", nothing) == "true"),
    authors = "Luiz M. Faria",
    sitename = "HierarchicalMatrices.jl",
    pages = Any["Home" => "index.md"
                "Modules" => "modules.md" ]
    # strict = true,
    # clean = true,
    # checkdocs = :exports,
)

deploydocs(
    repo = "github.com/maltezfaria/HierarchicalMatrices.jl.git",## Geometry
```@autodocs
Modules = [HierarchicalMatrices.Geometry]
Order   = [:module, :function, :type]
```
)
