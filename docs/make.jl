using Documenter
using CamiXon

makedocs(;
    modules=[CamiXon],
    authors="= <walra356@planet.nl> and contributors",
    repo="github.com/walra356/CamiXon.jl.git",
    sitename="CamiXon.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        size_threshold_warn = 250000,
        size_threshold = 300000,
        repolink="https://walra356.github.io/CamiXon.jl",
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
