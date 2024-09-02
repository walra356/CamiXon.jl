using Documenter
using CamiXon

makedocs(;
    modules=[CamiXon],
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
    ],
)

deploydocs(;
    repo="github.com/walra356/CamiXon.jl.git",
    devbranch = "main"
)