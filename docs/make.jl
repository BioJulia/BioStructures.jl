using Documenter, BioStructures

makedocs()
deploydocs(
    deps = Deps.pip("mkdocs", "pygments", "mkdocs-material"),
    repo = "github.com/BioJulia/BioStructures.jl.git",
    julia = "0.6",
    osname = "linux",
)
