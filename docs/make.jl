using PauliOperators
using Documenter

DocMeta.setdocmeta!(PauliOperators, :DocTestSetup, :(using PauliOperators); recursive=true)

makedocs(;
    modules=[PauliOperators],
    authors="Nick Mayhall",
    repo="https://github.com/nmayhall-vt/PauliOperators.jl/blob/{commit}{path}#{line}",
    sitename="PauliOperators.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://nmayhall-vt.github.io/PauliOperators.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Types" => "types.md",
        "Funcitons" => "functions.md",
    ],
)

deploydocs(;
    repo="github.com/nmayhall-vt/PauliOperators.jl",
    devbranch="main",
)
