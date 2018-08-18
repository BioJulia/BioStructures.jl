# BioStructures.jl

**Latest Release:**

[![Latest Release](https://img.shields.io/github/release/BioJulia/BioStructures.jl.svg)](https://github.com/BioJulia/BioStructures.jl/releases/latest)
[![BioStructures](http://pkg.julialang.org/badges/BioStructures_0.6.svg)](http://pkg.julialang.org/?pkg=BioStructures)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/BioJulia/BioStructures.jl/blob/master/LICENSE.md)
[![Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://biojulia.github.io/BioStructures.jl/stable)
![BioJulia maintainer: jgreener64](https://img.shields.io/badge/BioJulia%20Maintainer-jgreener64-orange.svg)

**Development status:**

[![Build Status](https://travis-ci.org/BioJulia/BioStructures.jl.svg?branch=master)](https://travis-ci.org/BioJulia/BioStructures.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/ltynlacyj689ei1u/branch/master?svg=true)](https://ci.appveyor.com/project/jgreener64/biostructures-jl/branch/master)
[![codecov.io](http://codecov.io/github/BioJulia/BioStructures.jl/coverage.svg?branch=master)](http://codecov.io/github/BioJulia/BioStructures.jl?branch=master)
[![Documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://biojulia.github.io/BioStructures.jl/latest)

## Description

BioStructures provides functionality to read, write and manipulate
macromolecular structures, in particular proteins.
[Protein Data Bank](https://www.rcsb.org/pdb/home/home.do) (PDB) and mmCIF
format files can be read in to a hierarchical data structure. Basic spatial
calculations and functions to access the PDB are also provided.

## Installation

Install BioStructures from the Julia REPL:

```julia
using Pkg
Pkg.add("BioStructures")
```

If you are interested in the cutting edge of the development, please check out
the master branch to try new features before release.

## Contributing and Questions

We appreciate contributions from users including reporting bugs, fixing issues,
improving performance and adding new features.
Please go to the [contributing section of the documentation](http://biojulia.github.io/BioStructures.jl/latest/contributing)
for more information.

If you have a question about
contributing or using this package, you are encouraged to use the
[Bio category of the Julia discourse
site](https://discourse.julialang.org/c/domain/bio).
