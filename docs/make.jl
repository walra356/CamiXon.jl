using CamiXon
using Documenter
using LaTeXStrings
using Pkg
Pkg.add("FITSIO")
using FITSIO 

makedocs(;
    modules=[CamiXon],
    authors="= <walra356@planet.nl> and contributors",
    repo="https://github.com/walra356/CamiXon.jl/blob/{commit}{path}#L{line}",
    sitename="CamiXon.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://walra356.github.io/CamiXon.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/walra356/CamiXon.jl",
)
