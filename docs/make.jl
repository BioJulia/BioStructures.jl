using Documenter
using BioStructures

makedocs(
    sitename = "BioStructures.jl",
    pages = [
        "Home" => "index.md",
        "Documentation" => "documentation.md",
        "Examples" => "examples.md",
        "API" => "api.md",
    ],
)

deploydocs(repo = "github.com/BioJulia/BioStructures.jl.git", push_preview = true)
