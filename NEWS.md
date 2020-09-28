BioStructures.jl release notes
==============================

## v0.11.3 - Sep 2020

* Gzip support is added for reading and writing mmCIF files via the `gzip` keyword argument.

## v0.11.2 - Sep 2020

* Add `get` method for `MMCIFDict`.

## v0.11.1 - Sep 2020

* Fix bug in reading mmCIF data values containing a comment character.

## v0.11.0 - Jun 2020

* The required versions of BioSequences.jl and BioAlignments.jl are updated to v2 of each, with support for earlier versions being dropped. `AminoAcidSequence` is hence renamed to `LongAminoAcidSeq`.
* `threeletter_to_aa`, a lookup table of amino acids, is re-exported from BioSymbols.

## v0.10.1 - May 2020

* Change compatibility bounds for new DataFrames.jl release.

## v0.10.0 - Apr 2020

* Change keyword argument names `pdb_dir` to `dir` and `file_format` to `format` for `downloadpdb`, `downloadentirepdb`, `updatelocalpdb`, `downloadallobsoletepdb` and `retrievepdb`.
* Remove `readpdb`, which has the same functionality as `read`.
* API reference section, more docstrings, links to related software and interactive Bio3DView.jl examples in documentation.

## v0.9.4 - Apr 2020

* Change compatibility bounds for new RecipesBase.jl and CodecZlib.jl releases.

## v0.9.3 - Mar 2020

* Change compatibility bounds for new RecipesBase.jl release.

## v0.9.2 - Feb 2020

* Improvements to performance throughout the package. Some functions are made up to 5 times faster.

## v0.9.1 - Jan 2020

* Fix documentation build.

## v0.9.0 - Jan 2020

* A reader and writer are added for the MMTF file format, building on top of [MMTF.jl](https://github.com/BioJulia/MMTF.jl). The interface is the same as for PDB and mmCIF files, with files either being read into the standard hierarchical structure or a `MMTFDict`. Gzipped files are supported. PDB, mmCIF and MMTF files can be interconverted.
* The `expand_disordered` flag is added to `collectatoms`, `collectresidues`, `countatoms`, `countresidues`, `coordarray`, `writepdb`, `writemmcif`, `writemmtf` and `DataFrame`. It determines whether disordered atoms and residues are expanded to include all entries. By default it is `false` except for the output functions, i.e. the last four above, where it is `true` by default.
* The `pdbextension` dictionary is changed to remove leading dots in the values.
* Improved file writing of empty elements.
* Examples are split off into a separate section in the documentation.
* A benchmark suite is added to track performance.

## v0.8.0 - Dec 2019

* Superimposition of structural elements is supported using the Kabsch algorithm. New functions are `superimpose!`, `Transformation`, `applytransform!` and `applytransform`.
* `rmsd` and `displacements` carry out superimposition by default, with the relevant keyword arguments available. Setting `superimpose` to `false` prevents this. `rmsdatoms` and `dispatoms` respectively determine which atoms to calculate the property for.
* The trivial `allselector`, which selects all atoms or residues, is added.
* The backbone oxygen `"O"` is added to `backboneatomnames`.
* Compatible bounds of package dependencies are added to Project.toml.

## v0.7.0 - Oct 2019

* `MetaGraph` from MetaGraphs.jl is extended to create graphs of contacting elements in a molecular structure, giving access to all the graph analysis tools in LightGraphs.jl.
* `DataFrame` from DataFrames.jl is extended to allow creation of data frames from lists of atoms or residues.
* `pairalign` from BioAlignments.jl is extended to produce pairwise alignments from structural elements.
* `AminoAcidSequence` now takes any element type and has the `gaps` keyword argument.
* Documentation example of interoperability with NearestNeighbors.jl.
* Parametric types used more extensively internally.

## v0.6.0 - Sep 2019

* `collectatoms`, `collectresidues`, `collectchains` and `collectmodels` no longer run `sort` before returning the final list. The user can run an explicit `sort` themselves if desired. This change makes the functions faster and allows preservation of the element order.
* Speed up residue iteration.
* Documentation improvements.

## v0.5.1 - Aug 2019

* Fix `MMCIFDict` to always contain a `Dict{String, Vector{String}}` rather than a `Dict{String, Union{String, Vector{String}}}`, which includes making the `"data_"` tag a `Vector{String}`.
* More functions documented and documentation bugfixes.

## v0.5.0 - Aug 2019

* The mmCIF reader now returns `Array{String,1}` for dictionary values even when there is only a single component. This improves consistency.
* Documentation expanded with references to Bio3DView.jl and an extra example.
* Replace REQUIRE with Project.toml.
* Bugfix when reading truncated MODEL line in a PDB file.

## v0.4.0 - Sep 2018

* The `ContactMap` and `DistanceMap` types are introduced along with their supertype `SpatialMap`. `contactmap` is removed. Plot recipes are defined for visualisation of `ContactMap`s and `DistanceMap`s. `showcontactmap` provides a quick way to view a `ContactMap` in the terminal.
* Bug fix on downloading MMTF files.

## v0.3.0 - Aug 2018

* Code is now compatible with Julia v0.7 and v1.0. Support for earlier Julia versions is dropped.
* `downloadpdb` can now be given a function as the first argument, in which case the function is run with the downloaded filepath(s) as an argument and the file(s) are deleted afterwards.
* Improved function docstrings.

## v0.2.0 - Mar 2018

* A reader and writer is added for the mmCIF format, which has been the standard PDB archive format since 2014. mmCIF files can either be read into a hierarchical structure object or directly in as a dictionary. PDB and mmCIF files can be interconverted.
* `chainid` now returns a `String` instead of a `Char`. This allows multi-character chain IDs. This also changes `chainids`, `chain` and `chains`. Chains can be accessed by string (e.g. `struc["A"]`), but can still be accessed by character for single chain IDs (e.g. `struc['A']`).
* `show` now returns a single line statement for objects across the module, in line with Julia conventions.

## v0.1.0 - Aug 2017

Transfer of existing code from [Bio.jl](https://github.com/BioJulia/Bio.jl). Compatible with Julia v0.6.

Features:
* Hierarchical data structure suitable for macromolecules, particularly proteins.
* Fast reader and writer for the [Protein Data Bank](https://www.rcsb.org/pdb/home/home.do) (PDB) file format.
* Selection and iteration of structural elements.
* Calculation of spatial properties such as distances and Ramachandran angles.
* Functions to access the PDB.
