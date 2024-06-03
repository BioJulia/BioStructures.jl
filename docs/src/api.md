# BioStructures API

```@meta
CurrentModule = BioStructures
```

On Julia 1.9 and later package extensions are used in order to reduce the number of dependencies:
- To use `LongAA`, call `using BioSequences`.
- To use `pairalign` or `Transformation` on structural elements, call `using BioSequences, BioAlignments`.
- To use `DataFrame`, call `using DataFrames`.
- To use `MetaGraph`, call `using Graphs, MetaGraphs`.

Exported names:
```@index
Order   = [:module, :type, :constant, :function, :macro]
```

Non-exported names:
- [`BioStructures.x`](@ref)
- [`BioStructures.x!`](@ref)
- [`BioStructures.y`](@ref)
- [`BioStructures.y!`](@ref)
- [`BioStructures.z`](@ref)
- [`BioStructures.z!`](@ref)

Docstrings:
```@autodocs
Modules = [BioStructures]
Private = false
Order   = [:module, :type, :constant, :function, :macro]
```
```@docs
x
x!
y
y!
z
z!
```
