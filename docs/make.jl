using Documenter
using GeometrySDF

# Setup for doctests in docstrings
DocMeta.setdocmeta!(GeometrySDF, :DocTestSetup, :(using GeometrySDF); recursive=true)

makedocs(;
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
    ),
    modules = [GeometrySDF],
    sitename = "GeometrySDF.jl",
    pages=[
        "Home" => "index.md",
        "API" => "api.md",
    ],
    doctest = true, # :fix
    warnonly = [:missing_docs],
)

deploydocs(
    repo = "github.com/KeitaNakamura/GeometrySDF.jl.git",
    devbranch = "main",
    push_preview = true,
)
