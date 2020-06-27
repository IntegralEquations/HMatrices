push!(LOAD_PATH,joinpath(@__DIR__, ".."))
using Documenter, HMatrices

makedocs(
    modules = [HMatrices],
    format = Documenter.HTML(; prettyurls = get(ENV, "CI", nothing) == "true"),
    authors = "Luiz M. Faria",
    sitename = "HMatrices.jl",
    pages = Any["Home" => "index.md"
                "Modules" => "modules.md" ]
    # strict = true,
    # clean = true,
    # checkdocs = :exports,
)

deploydocs(
    repo = "github.com/maltezfaria/HMatrices.jl.git",## Geometry
```@autodocs
Modules = [HMatrices.Geometry]
Order   = [:module, :function, :type]
```
)
