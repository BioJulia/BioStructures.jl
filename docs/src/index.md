# BioStructures.jl

[![Build Status](https://travis-ci.org/BioJulia/BioStructures.jl.svg?branch=master)](https://travis-ci.org/BioJulia/BioStructures.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/ltynlacyj689ei1u/branch/master?svg=true)](https://ci.appveyor.com/project/jgreener64/biostructures-jl/branch/master)
[![codecov.io](http://codecov.io/github/BioJulia/BioStructures.jl/coverage.svg?branch=master)](http://codecov.io/github/BioJulia/BioStructures.jl?branch=master)
[![Documentation](https://img.shields.io/badge/docs-latest-blue.svg)](http://biojulia.net/Bio.jl/latest/man/structure)

## Description

BioStructures provides functionality to read, write and manipulate
macromolecular structures, in particular proteins.
[Protein Data Bank](https://www.rcsb.org/pdb/home/home.do) (PDB) format files
can be read in to a heirarchical data structure. Basic spatial calculations and
functions to access the PDB are also provided.

## Installation

Install BioStructures from the Julia REPL:

```julia
julia> Pkg.clone("https://github.com/BioJulia/BioStructures.jl")
```

If you are interested in the cutting edge of the development, please check out
the master branch to try new features before release.
