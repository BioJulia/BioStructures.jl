# BioStructures.jl
# ================
#
# A Julia package to read, write and manipulate macromolecular structures.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/BioStructures.jl/blob/master/LICENSE.md

"Read, write and manipulate macromolecular structures."
module BioStructures

using LinearAlgebra
using Statistics

using BioAlignments
using BioGenerics
using BioSequences
using BioSymbols
using CodecZlib
using Downloads
using Format
import MMTF: parsemmtf, writemmtf # Imported to avoid clash with MMTF name
using RecipesBase
using STRIDE_jll
using DSSP_jll

include("model.jl")
include("pdb.jl")
include("mmcif.jl")
include("mmtf.jl")
include("spatial.jl")
include("secondary.jl")

if !isdefined(Base, :get_extension)
    include("../ext/BioStructuresDataFramesExt.jl")
    include("../ext/BioStructuresGraphsExt.jl")
end

using PrecompileTools

@compile_workload begin
    let
        mktemp() do path, io
            print(io, """
                    ATOM      1  N   MET A   1      -8.706  17.775  26.032  1.00 54.68           N
                    ATOM      2  CA  MET A   1      -7.260  17.734  25.709  1.00 54.68           C
                    ATOM      3  C   MET A   1      -6.961  18.942  24.830  1.00 54.68           C
                    ATOM      4  CB  MET A   1      -6.938  16.397  25.011  1.00 54.68           C
                    ATOM      5  O   MET A   1      -7.640  19.078  23.823  1.00 54.68           O
                    ATOM      6  CG  MET A   1      -5.440  16.103  24.920  1.00 54.68           C
                    ATOM      7  SD  MET A   1      -5.078  14.380  24.476  1.00 54.68           S
                    ATOM      8  CE  MET A   1      -3.296  14.391  24.798  1.00 54.68           C
                    END
                    """)
            read(path, PDB)
        end
    end
end

end # BioStructures
