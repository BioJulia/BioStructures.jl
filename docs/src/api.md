# BioStructures API

```@meta
CurrentModule = BioStructures
```

Package extensions are used in order to reduce the number of dependencies:
- To use `LongAA`, call `using BioSequences`.
- To use `pairalign`, [`superimpose!`](@ref), [`rmsd`](@ref)/[`displacements`](@ref) with the `superimpose` option or [`Transformation`](@ref) on structural elements, call `using BioSequences, BioAlignments`.
- To use `DataFrame`, call `using DataFrames`.
- To use `MetaGraph`, call `using Graphs, MetaGraphs`.
- To use [`MMTFDict`](@ref) or [`writemmtf`](@ref), call `import MMTF`.
- To use [`rundssp!`](@ref), [`rundssp`](@ref) or the `run_dssp` option for `read`/[`retrievepdb`](@ref), call `using DSSP_jll`.
- To use [`runstride!`](@ref), [`runstride`](@ref) or the `run_stride` option for `read`/[`retrievepdb`](@ref), call `using STRIDE_jll`.

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
