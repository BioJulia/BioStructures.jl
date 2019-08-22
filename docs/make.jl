using Documenter, BioStructures

makedocs(
    format = :html,
    sitename = "BioStructures.jl",
    pages = [
        "Home"         => "index.md",
        "Documentation"=> "documentation.md"
    ],
    authors = "Joe G Greener, The BioJulia Organisation and other contributors."
)

deploydocs(
    repo = "github.com/BioJulia/BioStructures.jl.git",
    julia = "1.0",
    osname = "linux",
    target = "build",
    deps = nothing,
    make = nothing
)
