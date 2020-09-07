# BioStructures.jl

**Latest Release:**

[![Latest release](https://img.shields.io/github/release/BioJulia/BioStructures.jl.svg)](https://github.com/BioJulia/BioStructures.jl/releases/latest)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/BioJulia/BioStructures.jl/blob/master/LICENSE.md)
[![Documentation stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://biojulia.github.io/BioStructures.jl/stable)
![BioJulia maintainer: jgreener64](https://img.shields.io/badge/BioJulia%20Maintainer-jgreener64-orange.svg)

**Development status:**

[![Build status](https://travis-ci.org/BioJulia/BioStructures.jl.svg?branch=master)](https://travis-ci.org/BioJulia/BioStructures.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/ltynlacyj689ei1u/branch/master?svg=true)](https://ci.appveyor.com/project/jgreener64/biostructures-jl/branch/master)
[![codecov.io](http://codecov.io/github/BioJulia/BioStructures.jl/coverage.svg?branch=master)](http://codecov.io/github/BioJulia/BioStructures.jl?branch=master)
[![Documentation dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://biojulia.github.io/BioStructures.jl/dev)

## Description

BioStructures provides functionality to read, write and manipulate
macromolecular structures, in particular proteins.
[Protein Data Bank](https://www.rcsb.org/pdb/home/home.do) (PDB), mmCIF and MMTF
format files can be read in to a hierarchical data structure. Spatial
calculations and functions to access the PDB are also provided.
It compares favourably in terms of performance to other PDB parsers -
see some [benchmarks online](https://github.com/jgreener64/pdb-benchmarks).

## Installation

Install BioStructures from the Julia package REPL, which can be accessed by
pressing `]` from the Julia REPL:

```
add BioStructures
```

See the [documentation](https://biojulia.github.io/BioStructures.jl/stable) for information on how
to use BioStructures.

## Citation

If you use BioStructures, please cite the paper:

Greener JG, Selvaraj J and Ward BJ. BioStructures.jl: read, write and manipulate macromolecular structures in Julia, *Bioinformatics* (2020) - [link](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/btaa502/5837108?guestAccessKey=aec90643-1d43-4521-9883-4a4a669187da)

## Contributing and questions

We appreciate contributions from users including reporting bugs, fixing issues,
improving performance and adding new features.

Detailed guidance for contributing to all BioJulia packages is provided at the
[BioJulia Contribution Documentation](https://github.com/BioJulia/BioStructures.jl/blob/master/CONTRIBUTING.md).

If you have a question about contributing or using this package, you are
encouraged to use the
[BioJulia Gitter](https://gitter.im/BioJulia/General) or the
[Bio category of the Julia discourse site](https://discourse.julialang.org/c/domain/bio).
