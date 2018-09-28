BioStructures.jl release notes
==============================

## v0.4.0 - September 2018

* The `ContactMap` and `DistanceMap` types are introduced along with their supertype `SpatialMap`. `contactmap` is removed. Plot recipes are defined for visualisation of `ContactMap`s and `DistanceMap`s. `showcontactmap` provides a quick way to view a `ContactMap` in the terminal.
* Bug fix on downloading MMTF files.

## v0.3.0 - August 2018

* Code is now compatible with Julia v0.7 and v1.0. Support for earlier Julia versions is dropped.
* `downloadpdb` can now be given a function as the first argument, in which case the function is run with the downloaded filepath(s) as an argument and the file(s) are deleted afterwards.
* Improved function docstrings.

## v0.2.0 - March 2018

* A reader and writer is added for the mmCIF format, which has been the standard PDB archive format since 2014. mmCIF files can either be read into a hierarchical structure object or directly in as a dictionary. PDB and mmCIF files can be interconverted.
* `chainid` now returns a `String` instead of a `Char`. This allows multi-character chain IDs. This also changes `chainids`, `chain` and `chains`. Chains can be accessed by string (e.g. `struc["A"]`), but can still be accessed by character for single chain IDs (e.g. `struc['A']`).
* `show` now returns a single line statement for objects across the module, in line with Julia conventions.

## v0.1.0 - August 2017

Transfer of existing code from [Bio.jl](https://github.com/BioJulia/Bio.jl). Compatible with Julia v0.6.

Features:
* Hierarchical data structure suitable for macromolecules, particularly proteins.
* Fast reader and writer for the [Protein Data Bank](https://www.rcsb.org/pdb/home/home.do) (PDB) file format.
* Selection and iteration of structural elements.
* Calculation of spatial properties such as distances and Ramachandran angles.
* Functions to access the PDB.
