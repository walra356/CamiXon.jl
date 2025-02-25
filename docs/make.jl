using Documenter
using DocumenterInterLinks
using CamiXon
using CamiDiff

links = InterLinks(
    "CamiMath" => "https://walra356.github.io/CamiMath.jl/stable/",
    "CamiDiff" => "https://walra356.github.io/CamiDiff.jl/stable/"
)

makedocs(;
    modules=[CamiXon],
    plugins=[links],
    authors="= <walra356@planet.nl> and contributors",
    sitename="CamiXon.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        size_threshold_warn = 250000,
        size_threshold = 300000,
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Atom" => "man/properties/atom.md",
        "Coulomb integrals" => "man/coulomb.md",
        "Index" => "man/index.md"
    ],
)

deploydocs(;
    repo="github.com/walra356/CamiXon.jl.git",
    devbranch = "main"
)