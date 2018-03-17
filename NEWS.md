BioStructures.jl v0.2.0 Release Notes
=====================================

* A reader and writer is added for the mmCIF format, which has been the standard PDB archive format since 2014. mmCIF files can either be read into a hierarchical structure object or directly in as a dictionary. PDB and mmCIF files can be interconverted.
* `chainid` now returns a `String` instead of a `Char`. This allows multi-character chain IDs. This also changes `chainids`, `chain` and `chains`. Chains can be accessed by string (e.g. `struc["A"]`), but can still be accessed by character for single chain IDs (e.g. `struc['A']`).
* `show` now returns a single line statement for objects across the module, in line with Julia conventions.


BioStructures.jl v0.1.0 Release Notes
=====================================

Transfer of existing code from [Bio.jl](https://github.com/BioJulia/Bio.jl). Compatible with Julia v0.6.

Features
--------

* Hierarchical data structure suitable for macromolecules, particularly proteins.
* Fast reader and writer for the [Protein Data Bank](https://www.rcsb.org/pdb/home/home.do) (PDB) file format.
* Selection and iteration of structural elements.
* Calculation of spatial properties such as distances and Ramachandran angles.
* Functions to access the PDB.
