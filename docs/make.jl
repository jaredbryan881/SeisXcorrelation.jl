push!(LOAD_PATH,"../src/")

include("../src/SeisXcorrelation.jl")

using Documenter, .SeisXcorrelation

makedocs(;
    modules=[SeisXcorrelation],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/kura-okubo/SeisXcorrelation.jl/blob/{commit}{path}#L{line}",
    sitename="SeisXcorrelation.jl",
    authors="kurama",
)

deploydocs(;
    repo="github.com/kura-okubo/SeisXcorrelation.jl",
    target = "build",
    deps   = nothing,
    make   = nothing,
)
